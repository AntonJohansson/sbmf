#pragma once

#include <sbmf/types.h>
#include <pthread.h>

/*
 * Functions that are to be run as work
 * in the thread pool are required to
 * fulfill the following prototype.
 */
typedef void work_func(void*);

struct work {
	work_func* func;
	void* arg;
};

enum thread_pool_priority {
		LOW_PRIO = 0,
		HIGH_PRIO
};

struct thread_pool;

/*
 * Creates a new pool with "thread_count" threads. Threads are
 * started immediately and ready to accept work.
 */
struct thread_pool* thread_pool_create(u8 thread_count);

/*
 * Adds function + arguments to be executed by a thread whenever one
 * is available.
 */
void thread_pool_push_work(struct thread_pool* pool,
		enum thread_pool_priority priority,
		work_func* work_func, void* arg);

void thread_pool_push_work_stack(struct thread_pool* pool,
		enum thread_pool_priority priority,
		struct work* w, size_t size);

/*
 * Destroys the pool itself, thus invalidating the thread_pool* passed in.
 * Freeing of the pool will only happen once all threads are idle and
 * there is no more work to be done.
 */
void thread_pool_destroy(struct thread_pool* pool);

/*
 * Returns number the size of the work queue, i.e.
 * the amount of work backlogged.
 */
u64 thread_pool_queue_size(struct thread_pool* pool, enum thread_pool_priority priority);

void thread_pool_caller_do_work(struct thread_pool* pool);
void thread_pool_auto_signal(struct thread_pool* pool, bool b);
void thread_pool_signal(struct thread_pool* pool);
void thread_pool_broadcast(struct thread_pool* pool);
