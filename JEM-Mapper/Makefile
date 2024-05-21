CC = mpic++
CFLAGS = -std=c++14 -O3 -fopenmp -DLMER_LENGTH=8 
INCLUDE = -I$(PWD)/includes

GET_WINSIZE=1
KVAL=$(shell echo $(ksize)-$(GET_WINSIZE) | bc)

ifeq ($(shell expr $(KVAL) \<= 0), 1)
      KVAL=31
endif

ifeq ($(shell expr $(KVAL) \<= 31), 1)
  CFLAGS += -DWINDW_SIZE=$(KVAL)
endif



 



SRC = $(wildcard src/*.cpp)
HDR = $(wildcard includes/*.h)
NAME = jem

all: $(NAME)

$(NAME) : $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o $(NAME) $? -fopenmp

clean:
	rm -f $(NAME)

distclean: clean
