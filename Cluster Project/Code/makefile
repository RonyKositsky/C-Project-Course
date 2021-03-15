OBJS	= main.o algorithms.o group.o pIteration.o spmat.o eHandler.o
SOURCE	= main.c algorithms.c group.c pIteration.c spmat.c eHandler.c
HEADER	= algorithms.h group.h pIteration.h spmat.h eHandler.h
OUT	= cluster
CC	 = gcc
FLAGS	 = -g -c -ansi -Wall -Wextra -Werror -pedantic-errors
LFLAGS	 = -lm

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.c
	$(CC) $(FLAGS) main.c 

algorithms.o: algorithms.c
	$(CC) $(FLAGS) algorithms.c 

group.o: group.c
	$(CC) $(FLAGS) group.c 

pIteration.o: pIteration.c
	$(CC) $(FLAGS) pIteration.c 

spmat.o: spmat.c
	$(CC) $(FLAGS) spmat.c 

eHandler.o: eHandler.c
	$(CC) $(FLAGS) eHandler.c 


clean:
	rm -f $(OBJS) $(OUT)