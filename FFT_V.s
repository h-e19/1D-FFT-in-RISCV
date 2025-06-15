#define STDOUT 0xd0580000
.set vectorSize, 128      # Change this for max num of 32 bit vector element supported by hardware

.section .text
.global _start
_start:

main:     
    li a0, vectorSize/2
    la a1, T_real
    la a2, T_imag
    call twiddle_gen
    la a0, vectorSize
    la a1, signal_re
    la a2, signal_im
    call bitrev_and_reorder_V


    la a0, signal_re
    la a1, signal_im
    la a2, T_real
    la a3, T_imag
    li a4, vectorSize
    call fft_V 

    la a0, signal_re
    li a1, vectorSize
    call printToLogVectorized
    la a0, signal_im
    li a1, vectorSize
    call printToLogVectorized

    j _finish                       # End program
    
    
logInt:                     # Takes input N(a0), returns its log-base 2 in a0
    addi sp, sp, -4                 # Make space to save registers used
    sw t0, 0(sp)                    # Save t0 as we will be using it
    
    add t0, a0, zero                # t0 = k = N    , loop counter
    add a0, zero, zero              # a0 = i = 0

    logLoop:                        # while(k) i.e when k is non-zero
    beq t0, zero, logLoopEnd        # If k becomes 0, end loop
    srai t0, t0, 1                  # Shift right k by 1. K = k >> 1
    addi a0, a0, 1                  # i++. Increment i
    j logLoop
    logLoopEnd:

    addi a0, a0, -1                 # Return i - 1
    lw t0, 0(sp)                    # Restore t0
    addi sp, sp, 4                  # Restore stack pointer
    jr ra                           # Return to caller


power2:         # a0 = n
    li    a1, 1       
    sll   a1, a1, a0  
    jr ra

twiddle_gen: # a0 = N/2, a1 = T_re addr, a2 = T_im addr
    addi sp, sp, -16
    sw ra, 0(sp)
    sw a0, 4(sp)
    sw a1, 8(sp)
    sw a2, 12(sp)
#----------------------------------angle processing-----------------------------------------------
    lw a0, 4(sp)
    vsetvli t0, a0, e32, m1
    vle32.v v0, (a1)            #v0 = T_real
    vle32.v v1, (a2)            #v1 = T_imag
    
    la t1, TWO_PI
    flw ft0, 0(t1)              #ft0 = 2*pi
    mv t2, a0 
    slli t2, t2, 1                  
    fcvt.s.w  ft1, t2           #ft1 = N

    vfmul.vf v0, v0, ft0
    vfmul.vf v1, v1, ft0
    vfdiv.vf v0, v0, ft1
    vfdiv.vf v1, v1, ft1

    vse32.v v0, (a1)
    vse32.v v1, (a2)
#-----------------------------------------------------------------------------------------------

                                #v0 = x
    vfmul.vv v2, v0, v0         #v2 = x^2
    vfmul.vv v3, v2, v0         #v3 = x^3
    vfmul.vv v4, v3, v0         #v4 = x^4
    vfmul.vv v5, v4, v0         #v5 = x^5
    vfmul.vv v6, v5, v0         #v6 = x^6
    vfmul.vv v7, v6, v0         #v7 = x^7


#-------------------------------sin aprox-------------------------------------------------------

    la t0, f_sin_2nd
    flw ft0, 0(t0)
    vfmul.vf v3, v3, ft0        #v3 = (x^3)*1/3!

    la t0, f_sin_3rd
    flw ft0, 0(t0)
    vfmul.vf v5, v5, ft0        #v5 = (x^5)*1/5!

    la t0, f_sin_4th
    flw ft0, 0(t0)              
    vfmul.vf v7, v7, ft0        #v7 = (x^7)*1/7!

    vfsub.vv v0, v3, v0
    vfsub.vv v0, v0, v5
    vfadd.vv v0, v0, v7

    vse32.v v0, (a2)
#----------------------------------------------------------------------------------------------
#-------------------------------cos aprox------------------------------------------------------

    la t0, f_cos_2nd
    flw ft0, 0(t0)
    vfmul.vf v2, v2, ft0        #v2 = (x^2)*-1/2!

    la t0, f_cos_3rd
    flw ft0, 0(t0)
    vfmul.vf v4, v4, ft0        #v4 = (x^4)*1/4!

    la t0, f_cos_4th
    flw ft0, 0(t0)              
    vfmul.vf v6, v6, ft0        #v6 = (x^6)*1/6!

    la t0, f_one
    flw ft0, 0(t0)

    vfadd.vf v1, v2, ft0
    vfadd.vv v1, v1, v4
    vfsub.vv v1, v1, v6

    vse32.v v1, (a1)

    lw ra, 0(sp)
    addi sp, sp, 16
    ret

bitrev_and_reorder_V: # a0 = size, a1 = &x[]_re, a2 = &x[]_im

    addi sp, sp, -24
    sw ra, 0(sp)
    sw a0, 4(sp) #log N
    sw a1, 8(sp) # &x[]_re
    sw a2, 12(sp) # &x[]_im
    sw s0, 16(sp) # N
    sw s1, 20(sp) #shift size

    mv s0, a0 #s0 = size of array

    call logInt #a0 = logN
    sw a0, 8(sp)

    la t0, revtemp_real
    la t1, revtemp_imag
    
    #shift size according to logN
    li s1, 30
    sub s1, s1, a0 # total bits - our bit size = shiftsize
    
    vsetvli t2, s0, e32
    vid.v v1 #generate seq from 1 to vl-1

    li t3, 0 #num of elements done

    #load pattern bits for bit reversal
    li a3, 0x55555555
    li a4, 0x33333333
    li a5, 0x0F0F0F0F
    li a6, 0x00FF00FF

    bitrev_loop:
        bge t3, s0, end_bitrev_loop

        #bit rev for 32 bits - swap groups of 1, then 2, then 4, then 8, then 16

        #groups of 1 (single elements so odd and even)
        vsrl.vi v2, v1, 1   #shift right by 1
        vand.vx v2, v2, a3  #AND to mask (keep every other bit)
        vand.vx v3, v1, a3  #AND to mask (keep every other bit)
        vsll.vi v3, v3, 1   #shift left by 1
        vor.vv v2, v2, v3   #OR to store both even and odds back in same vector

        #groups of 2 
        vsrl.vi v3, v2, 2    #shift right by 2
        vand.vx v3, v3, a4  #AND to mask (keep every other pair)
        vand.vx v2, v2, a4  #AND to mask (keep every other bit)
        vsll.vi v2, v2, 2   #shift left by 2
        vor.vv v2, v2, v3   #OR to store all pairs in same vector

        #groups of 4 
        vsrl.vi v3, v2, 4    #shift right by 4
        vand.vx v3, v3, a5  #AND to mask (keep every other grp of 4)
        vand.vx v2, v2, a5  
        vsll.vi v2, v2, 4   #shift left by 4
        vor.vv v2, v2, v3   #OR to store all grps in same vector

        #groups of 8 
        vsrl.vi v3, v2, 8    #shift right by 8
        vand.vx v3, v3, a6  #AND to mask (keep every other grp of 8)
        vand.vx v2, v2, a6  
        vsll.vi v2, v2, 8   #shift left by 8
        vor.vv v2, v2, v3   #OR to store all grps in same vector

        #groups of 16 
        vsrl.vi v3, v2, 16    #shift right by 16
        vsll.vi v2, v2, 16   #shift left by 16
        vor.vv v2, v2, v3   #OR to store all grps in same vector


        #shift by bit size s1
        vsrl.vx v2, v2, s1

        #store reordered values & inc index in temp array
        vloxei32.v v4, 0(a1), v2
        vse32.v v4, 0(t0) 
        vloxei32.v v4, 0(a2), v2     
        vse32.v v4, 0(t1)

        slli t4, t2, 2 # 4 bytes per element
        add t0, t0, t4 
        add t1, t1, t4 

        #increment all values in series in v1 by vl
        vadd.vx v1, v1, t2

        add t3, t3, t2 #index of num of elements done 
        j bitrev_loop

    end_bitrev_loop:

    la t0, revtemp_real
    la t1, revtemp_imag
    li t3, 0

    reorder_loop:
        bge t3, s0, reorder_loop_end

        #load
        vle32.v v5, (t0) #revtemp_real
        vle32.v v6, (t1) #revtemp_imag

        #store
        vse32.v v5, (a1)
        vse32.v v6, (a2)

        slli t4, t2, 2
        add t0, t0, t4
        add t1, t1, t4
        add a1, a1, t4
        add a2, a2, t4
        add t3, t3, t2

        j reorder_loop
    reorder_loop_end:

    lw ra, 0(sp)
    lw a0, 4(sp) #log N
    lw a1, 8(sp) # &x[]_re
    lw a2, 12(sp) # &x[]_im
    lw s0, 16(sp) # N
    lw s1, 20(sp) #shift size
    addi sp, sp, 24

    jr ra

fXk_l: #a0 = &x[], a1 = k, returns x[k] in f0
    slli a1, a1, 2
    add a0, a0, a1
    flw f0, 0(a0)
    jr ra

fXk_s: #a0 = &x[], a1 = k, f0 = float to store
    slli a1, a1, 2
    add a0, a0, a1
    fsw f0, 0(a0)
    jr ra

fft_V:  #a0=&x[]_re, a1=&x[]_im, a2=&T[]_re, a3=&T[]_im, a4=N or size
    addi sp, sp, 24
    sw a0, 0(sp)
    sw a1, 4(sp)
    sw a2, 8(sp)
    sw a3, 12(sp)
    sw a4, 16(sp)
    sw ra, 20(sp)

    mv a0, a4
    call logInt
    mv s0, a0               #s0 = logN
    li s1, 1                #s1 = stage

    fft_V_stage:
        bgt s1, s0, fft_V_stage_end
        li s2, 1
        sll s2, s2, s1      #s2 = butterfly_size
        srl s3, s2, 1       #s3 = butterfly_split
        lw s4, 16(sp)       #s4 = N
        srl s4, s4, s1      #s4 = num_of_butterflies
        li s5, 1            #s5 = butterfly_num
        li s10, 0
        slli s9, s2, 2      #constant butterfly_size * 4 offset
        li a4, 0

        #s6 = Twiddle_stride; stride = 2 ^ (logN - stage)
        sub s7, s0, s1
        li s6, 1
        sll s6, s6, s7              #s6 = twiddle_stride
        slli s6, s6, 2

        vsetvli t0, s3, e32, m1

        lw t5, 8(sp)                #t5 = &T[]_re
        lw t6, 12(sp)               #t6 = &T[]_im

        vlse32.v v5, (t5), s6         #v5 = T[]_re
        vlse32.v v6, (t6), s6         #v6 = T[]_im

        fft_V_butterfly_num:
            bgt s5, s4, fft_V_butterfly_num_end

            lw s11, 0(sp)                
            add s11, s11, s10           #s11 = base address for upper_re
            mv t1, s3
            slli t1, t1, 2              #t1 = offset for lower
            add t2, t1, s11             #t2 = base address for lower_re
            lw t3, 4(sp)                
            add t3, t3, s10             #t3 = base address for upper_im
            add t4, t1, t3              #t4 = base address for lower_im


            vle32.v v1, (s11)           #v1 = X[upper]_re
            vle32.v v2, (t2)            #v2 = X[lower]_re
            vle32.v v3, (t3)            #v3 = X[upper]_im
            vle32.v v4, (t4)            #v4 = X[lower]_im


            vfmul.vv v7, v5, v2         #v7 = T[]_re * X[lower]_re
            vfmul.vv v8, v6, v4         #v8 = T[]_im * X[lower]_im
            vfsub.vv v9, v7, v8         #v9 = T[]_re * X[lower]_re - T[]_im * X[lower]_im

            vfmul.vv v10, v5, v4        #v10 = T[]_re * X[lower]_im
            vfmul.vv v11, v6, v2        #v11 = T[]_im * X[lower]_re
            vfadd.vv v12, v10, v11      #v12 = T[]_re * X[lower]_im + T[]_im * X[lower]_re

            vfadd.vv v13, v1, v9        #v13 = X[upper]_re + (T[]_re * X[lower]_re - T[]_im * X[lower]_im)
            vfadd.vv v14, v3, v12       #v14 = X[upper]_im + (T[]_re * X[lower]_im + T[]_im * X[lower]_re)

            vse32.v v13, (s11)
            vse32.v v14, (t3)

            vfsub.vv v15, v1, v9        #v15 = X[upper]_re - (T[]_re * X[lower]_re - T[]_im * X[lower]_im)
            vfsub.vv v16, v3, v12       #v16 = X[upper]_im - (T[]_re * X[lower]_im + T[]_im * X[lower]_re)

            vse32.v v15, (t2)
            vse32.v v16, (t4)

            #butterfly number offset
            add s10, s10, s9

            addi s5, s5, 1
            j fft_V_butterfly_num
        fft_V_butterfly_num_end:
    
        addi s1, s1, 1 
        j fft_V_stage
    fft_V_stage_end:

    lw ra, 20(sp)
    addi sp, sp, -24
    jr ra

# Logs values from array in a0 into registers v1 for debugging and output.
# Inputs:
#   - a0: Base address of array
#   - a1: Size of array i.e. number of elements to log
printToLogVectorized:        
    addi sp, sp, -8
    sw a0, 0(sp)
    sw a1, 4(sp)

    loop:

        li t0, 0x123                 # Pattern for help in python script
        li t0, 0x456                 # Pattern for help in python script

        mv a1, a1                   # moving size to get it from log 
        li t0,0 		                # load i = 0
        vsetvli t3, a1, e32           # Set VLEN based on a1

        vle32.v v1, (a0)              # Load real[i] into v1

        li t0, 0x123                    # Pattern for help in python script
        li t0, 0x456                    # Pattern for help in python script
        
        slli t4, t3, 2            
        add a0, a0, t4           
        sub a1, a1, t3           
        
        # Loop if remaining elements > 0
        bne a1,x0, loop

    lw a0, 0(sp)
    lw a1, 4(sp)
    addi sp, sp, 8

	jr ra


_finish:
    li x3, 0xd0580000
    addi x5, x0, 0xff
    sb x5, 0(x3)
    beq x0, x0, _finish
.rept 100
    nop
.endr


.data  
#signal_re: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0
#signal_im: .float 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

#T_real: .float 1, 1, 1, 1, 1, 1, 1, 1
#T_imag: .float 1, 1, 1, 1, 1, 1, 1, 1

signal_re: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0

signal_im: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0

revtemp_real: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0

revtemp_imag: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0

T_real: .float 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0

T_imag: .float 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0

f_one:     .float 1.0
f_zero:    .float 0.0
f_minus1:   .float -1.0

f_sin_2nd:  .float 0.1666666666666667
f_sin_3rd:  .float 0.0083333333333333
f_sin_4th:  .float 1.984126984126984e-4

f_cos_2nd:  .float -0.5
f_cos_3rd:  .float 0.0416666666666667
f_cos_4th:  .float 0.0013888888888889

PI: .float 3.14159265358979323846
NEG_PI: .float -3.14159265358979323846
TWO_PI: .float 6.28318530717958647692
NEG_TWO_PI: .float -6.28318530717958647692
HALF_PI: .float 1.57079632679489661923
NEG_HALF_PI: .float -1.57079632679489661923
ONE: .float 1
TERMS: .word 12

half_pi_hi:    .float 1.57079637e+0  # Ï€/2 high part
half_pi_lo:    .float -4.37113883e-8 # Ï€/2 low part
const_2_pi:    .float 6.36619747e-1  # 2/Ï€
const_12582912: .float 12582912.0    # 1.5 * 2^23
cos_coeff_0:   .float 2.44677067e-5  # Coefficient for cosine
cos_coeff_1:   .float -1.38877297e-3
cos_coeff_2:   .float 4.16666567e-2
cos_coeff_3:   .float -5.00000000e-1
cos_coeff_4:   .float 1.00000000e+0
sin_coeff_0:   .float 2.86567956e-6  # Coefficient for sine
sin_coeff_1:   .float -1.98559923e-4
sin_coeff_2:   .float 8.33338592e-3
sin_coeff_3:   .float -1.66666672e-1

PI_OVER_4: .float 0.785398  # ≈ π/4
