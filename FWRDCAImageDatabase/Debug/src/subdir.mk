################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/FWRDCAImageDatabase.cpp 

OBJS += \
./src/FWRDCAImageDatabase.o 

CPP_DEPS += \
./src/FWRDCAImageDatabase.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I"/home/srmq/cppworkspace/FWRDCA/src" -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -lm -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


