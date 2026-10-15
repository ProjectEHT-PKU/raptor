# GNU make's built-in CC=cc must not override the HDF5 wrapper default.
ifeq ($(origin CC),default)
CC = $(if $(HDF5_CC),$(HDF5_CC),h5cc)
endif
CPPFLAGS += -I$(RAPTOR)/src -I$(CURDIR)
CFLAGS = -fopenmp  -std=c99  -O2 -lm -lgsl -Wall -Wno-unused-but-set-variable -MMD -MP
LDFLAGS = -fopenmp -lm -lgsl -lgslcblas

VPATH=$(RAPTOR)/src

OBJDIR=build

TARGET=RAPTOR

SOURCES=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c pol_emission.c tetrad.c model.c constants.c camera.c axion.c
OBJECTS := $(patsubst %.c,$(OBJDIR)/%.o,$(SOURCES))

.PHONY: all clean
all: $(TARGET)

-include $(OBJECTS:.o=.d)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJDIR):
	mkdir -p $@

$(OBJECTS): $(OBJDIR)/%.o: %.c | $(OBJDIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


clean:
	rm -f $(OBJECTS) $(OBJECTS:.o=.d) $(TARGET)
