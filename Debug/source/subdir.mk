################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/hushdebug.cpp \
../source/matarray.cpp \
../source/mpo.cpp \
../source/mps.cpp \
../source/qoperator.cpp \
../source/qstate.cpp 

OBJS += \
./source/hushdebug.o \
./source/matarray.o \
./source/mpo.o \
./source/mps.o \
./source/qoperator.o \
./source/qstate.o 

CPP_DEPS += \
./source/hushdebug.d \
./source/matarray.d \
./source/mpo.d \
./source/mps.d \
./source/qoperator.d \
./source/qstate.d 


# Each subdirectory must supply rules for building sources it contributes
source/%.o: ../source/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


