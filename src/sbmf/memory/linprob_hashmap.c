typedef u64 linprob_hashmap_hash_func(void* key);

/* djb2 by Dan Bernstein */
static u64 hash_djb2(void* key) {
    u64 hash = 5381;

	for (u32 i = 0; i < sizeof(key); ++i) {
        hash = ((hash << 5) + hash) + *((u8*)key); /* hash * 33 + c */
	}

    return hash;
}


struct linprob_hashmap {
	u64 size;
	u64 occupied_size;

	u32 key_size;
	u32 value_size;
	u32 pair_size;


	u8* empty_key;
	u8* memory;

	linprob_hashmap_hash_func* hash;
};





static inline u8* _pair(struct linprob_hashmap* map, u32 i) {
	return map->memory + i*map->pair_size;
}

static inline u8* _pair_key(struct linprob_hashmap* map, u32 i) {
	return _pair(map, i);
}

static inline u8* _pair_value(struct linprob_hashmap* map, u32 i) {
	return _pair(map, i) + map->key_size;
}






static inline struct linprob_hashmap* linprob_hashmap_new(u64 size, u32 key_size, u32 value_size, void* empty_key, linprob_hashmap_hash_func* hash) {
	/* Ensure size power of 2 */
	u64 size_pow2 = 1;
	while (size_pow2 < size)
		size_pow2 <<= 1;

	void* mem = sbmf_stack_push(
			sizeof(struct linprob_hashmap) +
			key_size +
			(key_size + value_size) * size_pow2);

	struct linprob_hashmap* map = mem;
	map->size = size_pow2;
	map->occupied_size = 0;

	map->key_size = key_size;
	map->value_size = value_size;
	map->pair_size = key_size + value_size;
	map->empty_key = (u8*)(map + 1);

	memcpy(map->empty_key, empty_key, key_size);

	map->memory = map->empty_key + map->key_size;
	map->hash = hash;

	/* Make all keys the empty key */
	for (u32 i = 0; i < map->size; ++i) {
		memcpy(_pair_key(map, i), map->empty_key, map->key_size);
	}

	return map;
}

static inline u64 _key_to_index(struct linprob_hashmap* map, void* key) {
	// Note if the size is a power of two
	// 	index % size = index & (size - 1)
	const u64 mask = map->size - 1;
	return map->hash(key) & mask;
}

static inline u64 _probe_next(struct linprob_hashmap* map, u64 index) {
	// Note if the size is a power of two
	// 	index % size = index & (size - 1)
	const u64 mask = map->size - 1;
	return (index + 1) & mask;
}

static inline u64 _diff(struct linprob_hashmap* map, u64 a, u64 b) {
	const u64 mask = map->size - 1;
	// adding map->size in order to handle b > a
	return (map->size + (a-b)) & mask;
}





static inline void linprob_hashmap_insert(struct linprob_hashmap* map, void* key, void* value) {
	if (map->occupied_size >= map->size) {
		sbmf_log_error("linprob_hashmap_insert: error out of space");
		return;
	}
	for (u64 i = _key_to_index(map, key); ; i = _probe_next(map, i)) {
		u8* pair_key   = _pair_key(map, i);
		u8* pair_value = _pair_value(map, i);

		if (memcmp(pair_key, map->empty_key, map->key_size) == 0) {
			memcpy(pair_key, 	key, 	map->key_size);
			memcpy(pair_value, 	value, 	map->value_size);
			map->occupied_size++;
			break;
		} else if (memcmp(pair_key, key, map->key_size) == 0) {
			memcpy(pair_value, value, map->value_size);
			break;
		}
	}
}

static inline u8* _find(struct linprob_hashmap* map, void* key, u64* index) {
	if (memcmp(key, map->empty_key, map->key_size) == 0)
		sbmf_log_error("linprob_hashmap_find: empty key should not be used");

	u32 iter = 0;
	for(u64 i = _key_to_index(map, key); ; i = _probe_next(map, i)){
		u8* pair_key = _pair_key(map, i);
		if (memcmp(pair_key, key, map->key_size) == 0) {
			if (index)
				*index = i;
			return (void*)pair_key;
		}
		else if (memcmp(pair_key, map->empty_key, map->key_size) == 0) {
			return NULL;
		}

		if (iter >= map->size) {
			return NULL;
		}

		iter++;
	}
}

static inline void* linprob_hashmap_at(struct linprob_hashmap* map, void* key) {
	u8* pair = _find(map, key, NULL);
	if (pair == NULL)
		return NULL;

	/* value comes right after key */
	return (void*)(pair + map->key_size);
}

static inline bool linprob_map_erase(struct linprob_hashmap* map, void* key) {
	u64 index_to_remove;
	u8* pair = _find(map, key, &index_to_remove);
	if (pair == NULL)
		return false;

	// we have to handle items that were placed via collision hashes before we remove.
	for(u64 i = _probe_next(map, index_to_remove); ; i = _probe_next(map, i)){
		u8* pair = _pair(map, i);
		u8* pair_key = _pair_key(map, i);

		if (memcmp(pair_key, map->empty_key, map->key_size) == 0) {
			memcpy(pair_key, map->empty_key, map->key_size);
			map->size--;
			return true;
		} else {
			u64 ideal_index = _key_to_index(map, pair_key);
			if (_diff(map, index_to_remove, ideal_index) < _diff(map, i, ideal_index)) {
				// this branch is taken if the element we're looking at is not in its ideal position,
				// and the position of the element to remove happens to be close to that position.
				// This is done to avoid putting holes in the the table between items with colliding
				// hashes.
				memcpy(_pair(map, index_to_remove), pair, map->pair_size);
				index_to_remove = i;
			}
		}
	}

	return false;
}
