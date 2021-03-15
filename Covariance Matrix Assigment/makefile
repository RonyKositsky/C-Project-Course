NAME = tester
FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: $(NAME)
clean:
	rm -rf *.o $(NAME)

$(NAME): $(NAME).o
	gcc $(NAME).o -o $(NAME) $(LIBS)
$(NAME).o: $(NAME).c
	gcc $(FLAGS) -c $(NAME).c
