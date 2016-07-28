CCFLAGS	:= -O3 -std=c99

all: matpow


matpow: matpow.c common.o
	$(CC) $(CCFLAGS) $^ -o $@


matpow.o: matpow.c
	$(CC) $(CCFLAGS) -c $< -o $@


%.o: %.c %.h
	$(CC) $(CCFLAGS) -c $< -o $@


clean:
	@rm -f matpow
	@rm -f *.o


@PHONY: all clean
