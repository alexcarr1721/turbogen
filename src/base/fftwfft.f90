module fftwfftc
    use iso_fortran_env
    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.f03"

    contains

    subroutine fftw_fft_1d(f, inverse)
        complex(C_FLOAT_COMPLEX), intent(inout)  :: f(:)
        logical,    optional            :: inverse
        logical                         :: backward
        type(C_PTR)                     :: fwd_plan, bwd_plan 
        integer(C_INT)                  :: n

        if ( present(inverse) ) then
            backward = inverse 
        else
            backward = .false.
        end if

        n = size(f,dim=1)
        if ( backward ) then
            ! Take inverse FFT 

            ! Make a plan
            bwd_plan = fftwf_plan_dft_1d(n, f, f, FFTW_BACKWARD, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_dft(bwd_plan, f, f)
            ! Scale result
            f = f/real(n,real32)
            ! Destroy plan
            call fftwf_destroy_plan(bwd_plan)
        else
            ! Take forward FFT

            ! Make a plan
            fwd_plan = fftwf_plan_dft_1d(n, f, f, FFTW_FORWARD, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_dft(fwd_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(fwd_plan)
        end if
    end subroutine fftw_fft_1d

    subroutine fftw_fft_2d(f, inverse)
        complex(C_FLOAT_COMPLEX), intent(inout)  :: f(:,:)
        logical,    optional            :: inverse
        logical                         :: backward
        type(C_PTR)                     :: fwd_plan, bwd_plan 
        integer(C_INT)                  :: n, m

        if ( present(inverse) ) then
            backward = inverse 
        else
            backward = .false.
        end if

        n = size(f,dim=1)
        m = size(f,dim=2)
        if ( backward ) then
            ! Take inverse FFT 

            ! Make a plan
            bwd_plan = fftwf_plan_dft_2d(m, n, f, f, FFTW_BACKWARD, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_dft(bwd_plan, f, f)
            ! Scale result
            f = f/real(n*m,real32)
            ! Destroy plan
            call fftwf_destroy_plan(bwd_plan)
        else
            ! Take forward FFT

            ! Make a plan
            fwd_plan = fftwf_plan_dft_2d(m, n, f, f, FFTW_FORWARD, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_dft(fwd_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(fwd_plan)
        end if
    end subroutine fftw_fft_2d

    subroutine fftw_r2r_1d(f, type, inverse)
        real(C_FLOAT), intent(inout)    :: f(:)                                     ! 1D complex array of data 
        character(len=*)                :: type
        logical, optional               :: inverse
        type(C_PTR)                     :: r2r_plan
        integer(C_INT)                  :: n
        logical                         :: backward

        ! Compute forward or inverse transform?
        if ( present(inverse) ) then
            backward = inverse
        else
            backward = .false.
        end if

        ! Size of transform
        n = size(f,dim=1)

        ! Go to type of transform 
        if ( trim(type) .eq. "DCT1" ) then 
            goto 100
        else if ( trim(type) .eq. "DCT2" ) then
            goto 110
        else if ( trim(type) .eq. "DCT3" ) then
            goto 120
        else if ( trim(type) .eq. "DCT4" ) then
            goto 130
        else if ( trim(type) .eq. "DST1" ) then
            goto 200
        else if ( trim(type) .eq. "DST2" ) then
            goto 210
        else if ( trim(type) .eq. "DST3" ) then
            goto 220
        else if ( trim(type) .eq. "DST4" ) then
            goto 230
        else
            ! Default to DCT1 
            goto 100
        end if

        100 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        110 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        120 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        130 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        200 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        210 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        220 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        230 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_1d(n, f, f, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        999 continue

    end subroutine fftw_r2r_1d

    subroutine fftw_r2r_2d(f, type_x2, type_x3, inverse)
        real(C_FLOAT), intent(inout)    :: f(:,:)                               ! 1D complex array of data 
        character(len=*), intent(in)    :: type_x2, type_x3
        logical, optional               :: inverse
        type(C_PTR)                     :: r2r_plan
        integer(C_INT)                  :: n, m
        logical                         :: backward

        ! Compute forward or inverse transform?
        if ( present(inverse) ) then
            backward = inverse
        else
            backward = .false.
        end if

        ! Size of transform
        n = size(f,dim=2)
        m = size(f,dim=1)
        ! Go to type of transform 
        if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DCT1") ) then 
            goto 100
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 110
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 120
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 130
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 140
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 150
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 160
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DCT1") ) then
            goto 170
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DCT2") ) then 
            goto 200
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 210
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 220
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 230
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 240
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 250
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 260
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DCT2") ) then
            goto 270
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DCT3") ) then 
            goto 300
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 310
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 320
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 330
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 340
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 350
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 360
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DCT3") ) then
            goto 370
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DCT4") ) then 
            goto 400
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 410
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 420
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 430
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 440
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 450
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 460
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DCT4") ) then
            goto 470
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DST1") ) then 
            goto 500
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 510
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 520
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 530
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 540
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 550
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 560
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DST1") ) then
            goto 570
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DST2") ) then 
            goto 600
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 610
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 620
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 630
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 640
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 650
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 660
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DST2") ) then
            goto 670
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DST3") ) then 
            goto 700
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 710
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 720
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 730
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 740
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 750
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 760
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DST3") ) then
            goto 770
        else if ( (trim(type_x3) .eq. "DCT1") .and. (trim(type_x2) .eq. "DST4") ) then 
            goto 800
        else if ( (trim(type_x3) .eq. "DCT2") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 810
        else if ( (trim(type_x3) .eq. "DCT3") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 820
        else if ( (trim(type_x3) .eq. "DCT4") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 830
        else if ( (trim(type_x3) .eq. "DST1") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 840
        else if ( (trim(type_x3) .eq. "DST2") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 850
        else if ( (trim(type_x3) .eq. "DST3") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 860
        else if ( (trim(type_x3) .eq. "DST4") .and. (trim(type_x2) .eq. "DST4") ) then
            goto 870
        else
            ! Default to DCT1 and DCT1
            goto 100
        end if

        100 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*(real(m,real32) - 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        110 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) - 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        120 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) - 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        130 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) - 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        140 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*(real(m,real32) - 1.0) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        150 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) - 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        160 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) - 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        170 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) - 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 


        200 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        210 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        220 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        230 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        240 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*real(m,real32) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        250 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        260 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        270 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999

        300 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        310 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        320 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        330 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        340 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*real(m,real32) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        350 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        360 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        370 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999

        400 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        410 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        420 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        430 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        440 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*real(m,real32) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        450 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        460 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        470 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_REDFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999

        500 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*(real(m,real32) + 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        510 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) + 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        520 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) + 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        530 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) + 1.0) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        540 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*(real(m,real32) + 1.0) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        550 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) + 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        560 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) + 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        570 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*(real(m,real32) + 1.0) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT00, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 


        600 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        610 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        620 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        630 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        640 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*real(m,real32) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        650 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        660 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        670 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT10, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999

       700 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        710 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        720 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        730 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        740 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*real(m,real32) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        750 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        760 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        770 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT01, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999

        800 continue 
        ! DCT1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) - 1.0) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT00, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        810 continue 
        ! DCT2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT10, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        820 continue 
        ! DCT3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT01, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        830 continue 
        ! DCT4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) ) 
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_REDFT11, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        840 continue 
        ! DST1
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*(real(n,real32) + 1.0) * 2.0*real(m,real32) )
            ! Destroy plan 
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT00, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        850 continue 
        ! DST2
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT10, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        860 continue 
        ! DST3
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT01, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999 

        870 continue 
        ! DST4
        if ( backward ) then                                                        
            ! Compute backward transform

            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Scale data
            f = f/( 2.0*real(n,real32) * 2.0*real(m,real32) )
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        else                                                                        
            ! Compute forward transform
            
            ! Create plan
            r2r_plan = fftwf_plan_r2r_2d(n, m, f, f, FFTW_RODFT11, FFTW_RODFT11, FFTW_ESTIMATE)
            ! Execute transform
            call fftwf_execute_r2r(r2r_plan, f, f)
            ! Destroy plan
            call fftwf_destroy_plan(r2r_plan)
        end if
        goto 999

        999 continue

    end subroutine fftw_r2r_2d

    subroutine fftw_fft_r2r_2d(f, type_x2, type_x3, inverse)
        complex(real32), intent(inout)  :: f(:,:)
        character(len=*), intent(in)    :: type_x2, type_x3 
        logical, optional               :: inverse 
        ! Local 
        logical                         :: backward 
        integer(int32)                  :: i, j, k, n, m
        real(real32), allocatable       :: f1_real(:), f1_imag(:)
        real(real32), allocatable       :: f2_real(:), f2_imag(:)

        if ( present(inverse) ) then
            backward = inverse 
        else 
            backward = .false. 
        end if

        n = size(f,dim=1)
        m = size(f,dim=2)
        allocate ( f1_real(n), f1_imag(n), f2_real(m), f2_imag(m) )
        if ( type_x2 .eq. "FFT" ) then 
            do j = 1,m 
                call fftw_fft_1d(f(:,j), inverse=backward)
            end do 
        else 
            do j = 1,m 
                f1_real = real(f(:,j))
                f1_imag = imag(f(:,j))
                call fftw_r2r_1d(f1_real, type_x2, inverse=backward)
                call fftw_r2r_1d(f1_imag, type_x2, inverse=backward)
                f(:,j) = cmplx(f1_real,f1_imag)
            end do
        end if
        if ( type_x3 .eq. "FFT" ) then 
            do i = 1,n 
                call fftw_fft_1d(f(i,:), inverse=backward)
            end do
        else 
            do i = 1,n
                f2_real = real(f(i,:))
                f2_imag = imag(f(i,:))
                call fftw_r2r_1d(f2_real, type_x3, inverse=backward)
                call fftw_r2r_1d(f2_imag, type_x3, inverse=backward)
                f(i,:) = cmplx(f2_real,f2_imag)
            end do
        end if
        deallocate ( f1_real, f1_imag, f2_real, f2_imag )

    end subroutine fftw_fft_r2r_2d

    subroutine fftw_call_1d_for_2d(f, type1, type2, inverse)
        !***********************************************************************
        ! Purpose: 
        !           Perform triginometric transforms using intel mkl                                                       
        !                                                                        
        ! Inputs:    _________________________________________________________  
        !           |___Variable___|_______________Purpose____________________|  
        !           |      f       |  Complex input array                     |
        !           |    inverse   |  Logical arguement (if .true. take       | 
        !           |              |  inverse fft). Default is .false.        |
        !           |______________|__________________________________________|  
        !          
        ! Outputs:   _________________________________________________________
        !           |___Variable____|_______________Purpose___________________|  
        !           |      f        |  Complex potential array (in-place)     |
        !           |_______________|_________________________________________|  
        !                                                                     
        !***********************************************************************
        use mpi
        implicit none
        complex(real32), intent(inout)          :: f(:,:,:)
        character(len=*), optional              :: type1
        character(len=*), optional              :: type2
        logical, intent(in), optional           :: inverse
        ! Local variables ******************************************************
        integer(int32)                  :: i, j, k, n, m, l
        logical                         :: backward
        character(len=180)              :: t_type1, t_type2
        real(real32), allocatable       :: g1(:), g2(:), h1(:), h2(:)

        l = size(f,dim=3)
        m = size(f,dim=2)
        n = size(f,dim=1)
        allocate ( g1(n), g2(n), h1(m), h2(m) )

        if ( present(type1) ) then
            t_type1 = type1
        else
            t_type1 = "DCT1"
        end if

        if ( present(type2) ) then
            t_type2 = type2
        else
            t_type2 = "DCT1"
        end if

        if ( present(inverse) ) then
            backward = inverse
        else
            backward = .false.
        end if

        do k = 1,l 
            ! Compute transform in first dimension 
            do j = 1,m 
                ! Assign real and imaginary data to real arrays 
                do i = 1,n 
                    g1(i) = real(f(i,j,k))
                    g2(i) = imag(f(i,j,k))
                end do
                ! Take transforms of g1 and g2 (selection of transforms)
                call fftw_r2r_1d(g1, trim(t_type1), backward)
                call fftw_r2r_1d(g2, trim(t_type1), backward)
                ! Assign back to f
                do i = 1,n 
                    f(i,j,k) = cmplx(g1(i),g2(i))
                end do
            end do
            ! Compute transform in second dimension
            do i = 1,n 
                ! Assign real and imaginary data to real arrays 
                do j = 1,m 
                    h1(j) = real(f(i,j,k))
                    h2(j) = imag(f(i,j,k))
                end do
                ! Take transforms of h1 and h2 (selection of transforms)
                call fftw_r2r_1d(h1, trim(t_type2), backward)
                call fftw_r2r_1d(h2, trim(t_type2), backward)
                ! Assign back to f 
                do j = 1,m 
                    f(i,j,k) = cmplx(h1(j),h2(j))
                end do
            end do
        end do
      
    end subroutine fftw_call_1d_for_2d

end module fftwfftc