OBJ=cpml/cpml.o fd/fd.o getpar/getpar.o field/field.o record/record.o record/isis.o material/material.o material/mktable.o model/model.o model/util.o param/param.o point_source/waveform.o
CUOBJ=cpml/cu_cpml.o field/cu_field.o fd/cu_fd.o material/cu_material.o point_source/cu_point_source.o record/cu_record.o slide/cu_slide.o thread/cu_thread.o
CC=/usr/usc/cuda/5.5/bin/nvcc
PTHREAD=-lpthread
CFLAGS=-g -G
CFLAGS=-O3 -DUSE_RESTRICT -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -rdc=true
%.o: %.cpp 
	${CC} ${CFLAGS} -c -o $@ $< 
%.o: %.cu
	${CC} ${CFLAGS} -c -o $@ $< 
sh2d: ${OBJ} ${CUOBJ} main.o
	${CC} --ptxas-options=-v -m64 ${CFLAGS} ${PTHREAD} -o $@ $^ 
clean:
	\rm -f *.o */*.o
