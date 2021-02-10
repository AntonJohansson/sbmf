/*
 * This file implements helper functions for dealing with indices into multi-dimensional
 * triangular arrays. For instance mapping F: (i,j,k,l) -> I such that I is an index
 * into a linear array[N] and range(F) = [0,N), AND you can permute the arguments to in
 * whichever way and still get the same index.
 *
 * Implementations for the 2d and 4d case are provided, they all follow the same pattern tho.
 */

/*
 * Sorts 4 values in ascending order using
 * 5 comparisons, see sorting networks.
 */

__host__ __device__
static inline void sort4_cuda(u32* arr) {
#define SWAP(x,y) if (arr[x] > arr[y]) { int tmp = arr[x]; arr[x] = arr[y]; arr[y] = tmp; }
    SWAP(0,1);
    SWAP(2,3);
    SWAP(0,2);
    SWAP(1,3);
    SWAP(1,2);
#undef SWAP
}


/*
 * In general, a method to construct the map F: (i,j,k,l, ...) -> I is by defining maps f1,f2,... for each index and settings F: (i,j,k,l) -> f1(i) + f2(j) + f3(k) + f4(l) +  ..., where the maps f1,... can be defined recursively as
 *
 * 		f1: i -> i
 * 		f2: j -> sum{a=[0,j)} f1(a+1)
 * 		f3: k -> sum{a=[0,k)} f2(a+1)
 * 		f4: l -> sum{a=[0,l)} f3(a+1)
 * 		  .
 * 		  .
 * 		  .
 *
 * there is also some very nice graphical intuition as to why this works, but I won't get into that. The maps f1,f2,f3,f4 also have a closed form so we can avoid for loops, probably holds for larger f's, let wolfram alpha figure it out.
 */

/*
 * What follows is the implementation of the closed forms of f1,f2,f3
 */

__host__ __device__
static inline u32 subindex2_cuda(u32 j) {
    return (j*(j+1))/2;
}

__host__ __device__
static inline u32 subindex3_cuda(u32 k) {
    return (k*(k+1)*(k+2))/6;
}

__host__ __device__
static inline u32 subindex4_cuda(u32 l) {
    return (l*(l+1)*(l+2)*(l+3))/24;
}

/*
 * I assume the above pattern for the functions would continue and that's pretty cool.
 * Next comes the final index calculation, including a sort to ensure that it works for
 * differently ordered indices.
 */

__host__ __device__
static inline u32 index4_cuda(u32 i, u32 j, u32 k, u32 l) {
    u32 arr[4] = {i, j, k, l};
    sort4_cuda(arr);

    return arr[0] + subindex2_cuda(arr[1]) + subindex3_cuda(arr[2]) + subindex4_cuda(arr[3]);
}

/*
 * Lastly it's nice to include a wrapper to calculate the size of
 * these linear arrays
 */

__host__ __device__
static inline u64 size4_cuda(u64 N) {
	return index4_cuda(N-1, N-1, N-1, N-1) + 1;
}
