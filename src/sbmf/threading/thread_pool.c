#include "thread_pool.h"
#include <sbmf/sbmf.h>
#include <sbmf/debug/log.h>

#include <errno.h>
#include <unistd.h>
#include <sys/syscall.h>

#define TP_WORK_STACK_LEN 10

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
struct work_queue {
	size_t bottom, top, size;
	struct work stack[TP_WORK_STACK_LEN];
};

static inline struct work pop_work_queue(struct work_queue* queue){
	struct work w = {0};

	if (queue->size > 0) {
		queue->size--;
		w = queue->stack[queue->bottom];
		queue->bottom = (queue->bottom+1) % TP_WORK_STACK_LEN;
	}

	return w;
}

static inline bool push_work_queue(struct work_queue* queue, struct work w){
	if (queue->size < TP_WORK_STACK_LEN) {
		queue->stack[queue->top] = w;
		queue->top = (queue->top + 1) % TP_WORK_STACK_LEN;
		queue->size++;
		return true;
	}

	return false;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
struct thread_pool {
	struct work_queue low_prio, high_prio;

	// Locking
	pthread_mutex_t mutex;
	pthread_cond_t cond;

	bool running;
	bool auto_signal;

	u64 thread_count;
	pthread_t threads[];
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
static inline void print_cond_wait_res(int res){
	switch(res){
		case 0: 			  log_info("COND_WAIT OK\n");     break;
		case EPERM:           log_error("EPERM\n");           break;
		case ENOTRECOVERABLE: log_error("ENOTRECOVERABLE\n"); break;
		case ETIMEDOUT:       log_error("ETIMEDOUT\n");       break;
		case EINVAL:          log_error("EINVAL\n");          break;
		case EOWNERDEAD:      log_error("EOWNERDEAD\n");      break;
		default: 			  log_error("UNKNOWN ERROR\n");   break;
	};
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
static void* thread_worker_func(void* args){
	struct thread_pool* pool = (struct thread_pool*)args;

	struct work w = {0};
	i32 res = 0;

	while (true){
		pthread_mutex_lock(&pool->mutex);

		while(pool->running && pool->high_prio.size == 0 && pool->low_prio.size == 0)
			res = pthread_cond_wait(&pool->cond, &pool->mutex);

		//print_cond_wait_res(res);
		(void)res;

		if(!pool->running && pool->low_prio.size == 0 && pool->high_prio.size == 0){
			//log_info("pool stopped running, thread exiting\n");
			pthread_mutex_unlock(&pool->mutex);
			break;
		}

		if (pool->high_prio.size > 0)
			w = pop_work_queue(&pool->high_prio);
		else if (pool->low_prio.size > 0)
			w = pop_work_queue(&pool->low_prio);

		pthread_mutex_unlock(&pool->mutex);

		if(w.func != 0){
			//pid_t tid = syscall(SYS_gettid);
			//log_info("thread %d calling work func %p with arg %p\n", tid, w.func, w.arg);
			w.func(w.arg);
			w = (struct work){0,0};
		}
	}

	return NULL;
}

static void caller_worker_func(struct thread_pool* pool) {
	struct work w = {0};
	i32 res = 0;
	while (true) {
		pthread_mutex_lock(&pool->mutex);

		//print_cond_wait_res(res);
		(void)res;

		if(pool->low_prio.size == 0 && pool->high_prio.size == 0){
			pthread_mutex_unlock(&pool->mutex);
			break;
		}

		if (pool->high_prio.size > 0)
			w = pop_work_queue(&pool->high_prio);
		else if (pool->low_prio.size > 0)
			w = pop_work_queue(&pool->low_prio);

		pthread_mutex_unlock(&pool->mutex);

		if(w.func != 0){
			//pid_t tid = syscall(SYS_gettid);
			//log_info("thread %d calling work func %p with arg %p\n", tid, w.func, w.arg);
			w.func(w.arg);
			w = (struct work){0,0};
		}
	}
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
static struct work_queue* get_work_queue(struct thread_pool* pool, enum thread_pool_priority prio) {
	switch(prio){
		case HIGH_PRIO: return &pool->high_prio;
		case LOW_PRIO: return &pool->low_prio;
		default: return 0;
	};
}

struct thread_pool* thread_pool_create(u8 thread_count){
	struct thread_pool* pool = (struct thread_pool*) sbmf_stack_push(sizeof(struct thread_pool) + sizeof(pthread_t)*thread_count);

	pool->high_prio = (struct work_queue){.bottom=0,.top=0,.size=0};
	pool->low_prio  = (struct work_queue){.bottom=0,.top=0,.size=0};

	pool->thread_count = thread_count;
	pool->running = true;
	pool->auto_signal = true;

	pthread_mutex_init(&pool->mutex, 0);
	pthread_cond_init(&pool->cond, 0);

	for(u8 i = 0; i < thread_count; i++){
		pthread_create(&pool->threads[i], NULL, thread_worker_func, pool);
	}

	return pool;
}

void thread_pool_push_work(struct thread_pool* pool, enum thread_pool_priority priority, work_func* work_func, void* arg){
	struct work w = (struct work){work_func, arg};
	thread_pool_push_work_stack(pool, priority, &w, 1);
}

void thread_pool_push_work_stack(struct thread_pool* pool,
		enum thread_pool_priority priority, struct work* w, size_t size) {
	pthread_mutex_lock(&pool->mutex);
	struct work_queue* q = get_work_queue(pool, priority);
	for (size_t i = 0; i < size; i++) {
		if (!push_work_queue(q, w[i])) {
			log_error("thread_pool %p: %s work queue full, omitting work.\n", (void*)pool, (priority == HIGH_PRIO) ? "high prio." : "low prio.");
			break;
		}
	}
	pthread_mutex_unlock(&pool->mutex);

	if (pool->auto_signal)
		pthread_cond_signal(&pool->cond);
}

void thread_pool_destroy(struct thread_pool* pool) {
	pthread_mutex_lock(&pool->mutex);
	pool->running = false;
	pthread_mutex_unlock(&pool->mutex);
	pthread_cond_broadcast(&pool->cond);

	for(uint8_t i = 0; i < pool->thread_count; i++){
		pthread_join(pool->threads[i], NULL);
	}
}

size_t thread_pool_queue_size(struct thread_pool* pool, enum thread_pool_priority priority) {
	return (priority == LOW_PRIO) ? pool->low_prio.size : pool->high_prio.size;
}

void thread_pool_caller_do_work(struct thread_pool* pool) {
	caller_worker_func(pool);
}

void thread_pool_auto_signal(struct thread_pool* pool, bool b) {
	pool->auto_signal = b;
}

void thread_pool_signal(struct thread_pool* pool) {
	pthread_cond_signal(&pool->cond);
}
void thread_pool_broadcast(struct thread_pool* pool) {
	pthread_cond_broadcast(&pool->cond);
}
