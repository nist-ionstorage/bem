/* JA 080930: Memory tracker for automatic heap deallocation for FastLap */

#ifdef __cplusplus
extern "C" {
#endif

void mtinit();
void mtadd(void *ptr);
void mtclear();

#ifdef __cplusplus
}
#endif
