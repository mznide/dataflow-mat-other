CCFLAGS	:= -O3 -std=c99

all: matvecmulti


matvecmulti: matvecmulti.c common.o
	$(CC) $(CCFLAGS) $^ -o $@


matvecmulti.o: matvecmulti.c
	$(CC) $(CCFLAGS) -c $< -o $@


%.o: %.c %.h
	$(CC) $(CCFLAGS) -c $< -o $@


clean:
	@rm -f matvecmulti
	@rm -f *.o


@PHONY: all clean
