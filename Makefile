CC ?= gcc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -pedantic -I.
LDFLAGS ?= -lm

SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
DEPS := $(OBJS:.o=.d)
TARGET := dvm_sim

.PHONY: all run clean distclean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

-include $(DEPS)

%.o: %.c
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(DEPS)

distclean: clean
	rm -f $(TARGET)
