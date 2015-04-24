################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../fit_formulas.c \
../main.c \
../mymatrices.c 

OBJS += \
./fit_formulas.o \
./main.o \
./mymatrices.o 

C_DEPS += \
./fit_formulas.d \
./main.d \
./mymatrices.d 


# Each subdirectory must supply rules for building sources it contributes
fit_formulas.o: ../fit_formulas.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -lm -MMD -MP -MF"$(@:%.o=%.d)" -MT"fit_formulas.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


