#include <cstdint>

#ifdef __aarch64__
static int64_t get_utime(){
	uint64_t tsc;
	asm volatile ("mrs %0, cntvct_el0" : "=r" (tsc));
	return tsc;
}
static double tick2second(uint64_t tick){
	auto frequency = []{
		uint64_t frq;
		asm volatile ("mrs %0, cntfrq_el0" : "=r" (frq));
		return frq;
	};
	static double invfreq = 1.0 / frequency();
	return invfreq * (double)tick;
}
#else
#  ifdef __APPLE__
#  include <sys/time.h>
static int64_t get_utime(){
	timeval tv;
	gettimeofday(&tv, NULL);

	return tv.tv_usec*1000ll + tv.tv_sec*1000000000ll;
}
#  elif defined __linux__
#  include <time.h>
static int64_t get_utime(){
	timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);

	return ts.tv_nsec + ts.tv_sec*1000000000ll;
}
#  endif
static double tick2second(uint64_t tick){
	return 1.e-9 * (double)tick;
}
#endif
