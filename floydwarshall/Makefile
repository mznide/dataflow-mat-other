CCFLAGS	:= -O3 -std=c99

all: floydwarshall


floydwarshall: floydwarshall.c common.o
	$(CC) $(CCFLAGS) $^ -o $@


floydwarshall.o: floydwarshall.c
	$(CC) $(CCFLAGS) -c $< -o $@


%.o: %.c %.h
	$(CC) $(CCFLAGS) -c $< -o $@


clean:
	@rm -f floydwarshall
	@rm -f *.o


@PHONY: all clean
