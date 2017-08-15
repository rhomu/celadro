/** Warp reduce, from:
 * https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/
 */

template<class T>
__inline__ __device__
double warpReduceSum(double val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val += __shfl_down(val, offset);
  return val;
}

template<class T>
__inline__ __device__
int blockReduceSum(T val)
{
  static __shared__ double shared[warpSize];
  int lane = threadIdx.x % warpSize;
  int wid  = threadIdx.x / warpSize;

  val = warpReduceSum(val);     // Each warp performs partial reduction

  if (lane==0) shared[wid]=val; // Write reduced value to shared memory

  __syncthreads();              // Wait for all partial reductions

  //read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

  if (wid==0) val = warpReduceSum(val); //Final reduce within first warp

  return val;
}
