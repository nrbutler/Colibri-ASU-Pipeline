CC = gcc

CFLAGS = -g -Wall -lcfitsio -lm

all: biasreduce darkreduce0 darkreduce flatreduce0 flatreduce imarith immed imreduce0 imreduce1 imreduce weight_clip calc_psf calc_psf_bg calc_psf_shift calc_phot psf_deconvolve psf_subtract

biasreduce: biasreduce.c
	$(CC) -o biasreduce biasreduce.c $(CFLAGS)

darkreduce0: darkreduce0.c
	$(CC) -o darkreduce0 darkreduce0.c $(CFLAGS)

darkreduce: darkreduce.c
	$(CC) -o darkreduce darkreduce.c $(CFLAGS)

flatreduce0: flatreduce0.c
	$(CC) -o flatreduce0 flatreduce0.c $(CFLAGS)

flatreduce: flatreduce.c
	$(CC) -o flatreduce flatreduce.c $(CFLAGS)

imarith: imarith.c
	$(CC) -o imarith imarith.c $(CFLAGS)

immed: immed.c
	$(CC) -o immed immed.c $(CFLAGS)

imreduce0: imreduce0.c
	$(CC) -o imreduce0 imreduce0.c $(CFLAGS)

imreduce1: imreduce1.c
	$(CC) -o imreduce1 imreduce1.c $(CFLAGS)

imreduce: imreduce.c
	$(CC) -o imreduce imreduce.c $(CFLAGS)

weight_clip: weight_clip.c
	$(CC) -o weight_clip weight_clip.c $(CFLAGS)

calc_psf: calc_psf.c
	$(CC) -o calc_psf calc_psf.c $(CFLAGS)

calc_psf_bg: calc_psf_bg.c
	$(CC) -o calc_psf_bg calc_psf_bg.c $(CFLAGS)

calc_psf_shift: calc_psf_shift.c
	$(CC) -o calc_psf_shift calc_psf_shift.c $(CFLAGS)

calc_phot: calc_phot.c
	$(CC) -o calc_phot calc_phot.c $(CFLAGS)

psf_deconvolve: psf_deconvolve.c
	$(CC) -o psf_deconvolve psf_deconvolve.c $(CFLAGS)

psf_subtract: psf_subtract.c
	$(CC) -o psf_subtract psf_subtract.c $(CFLAGS)

install:
	cp biasreduce darkreduce0 darkreduce flatreduce0 flatreduce imarith immed imreduce0 imreduce1 imreduce weight_clip calc_psf calc_psf_bg calc_psf_shift calc_phot psf_deconvolve psf_subtract ..

clean:
	$(RM) biasreduce darkreduce0 darkreduce flatreduce0 flatreduce imarith immed imreduce0 imreduce1 imreduce weight_clip calc_psf calc_psf_bg calc_psf_shift calc_phot psf_deconvolve psf_subtract
