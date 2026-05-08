! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module MLT

implicit none

private
public :: calc_MLT

contains

   !> Calculates the outputs of convective mixing length theory.
   !!
   !! @param MLT_option A string specifying which MLT option to use. Currently supported are Cox, Henyey, ML1, ML2, Mihalas. Note that 'TDC' is also a valid input and will return the Cox result. This is for use when falling back from TDC -> MLT, as Cox is the most-similar prescription to TDC.
   !! @param mixing_length_alpha The mixing length parameter.
   !! @param Henyey_MLT_nu_param The nu parameter in Henyey's MLT prescription.
   !! @param Henyey_MLT_y_param The y parameter in Henyey's MLT prescription.
   !! @param chiT dlnP/dlnT|rho
   !! @param chiRho dlnP/dlnRho|T
   !! @param Cp Specific heat at constant pressure (erg/g/K).
   !! @param grav The acceleration due to gravity (cm/s^2).
   !! @param Lambda The mixing length (cm).
   !! @param rho density (g/cm^3).
   !! @param T temperature (K).
   !! @param opacity opacity (cm^2/g)
   !! @param gradr The radiative temperature gradient dlnT/dlnP_{rad}
   !! @param grada The adiabatic temperature gradient dlnT/dlnP|s
   !! @param gradL The Ledoux temperature gradient dlnT/dlnP
   !! @param Gamma The convective efficiency parameter (output).
   !! @param gradT The temperature gradient dlnT/dlnP (output).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param conv_vel The convection speed (cm/s).
   !! @param D The chemical diffusion coefficient (cm^2/s).
   !! @param mixing_type Set to convective if convection operates (output).
   !! @param ierr Tracks errors (output).
   subroutine calc_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, &
                     Henyey_MLT_y_param, chiT, chiRho, Cp, grav, Lambda, rho, &
                     P, T, opacity, gradr, grada, gradL, Gamma, gradT, Y_face, &
                     conv_vel, D, mixing_type, max_conv_vel, ierr)
      use const_def, only: dp, clight, convective_mixing, crad, two_13, four_13, one_third
      use num_lib
      use utils_lib
      use auto_diff
      type(auto_diff_real_star_order1), intent(in) :: chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=*), intent(in) :: MLT_option
      real(dp), intent(in) :: mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, max_conv_vel
      type(auto_diff_real_star_order1), intent(out) :: Gamma, gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      real(dp) :: ff1, ff2, ff3, ff4
      type(auto_diff_real_star_order1) :: &
         Q, omega, a0, ff4_omega2_plus_1, A_1, A_2, &
         A_numerator, A_denom, A, Bcubed, delta, Zeta, &
         f, f0, f1, f2, radiative_conductivity, convective_conductivity, csound, &
         chi


      real(dp):: gamma_val
      
      real(dp) :: fn_val

      real(dp) :: Sigma, dSigma_dGamma
      real(dp) :: S
      real(dp) :: OmegaCM, dOmegaCM_dGamma
      real(dp) :: OmegaCGM, dOmegaCGM_dGamma
      real(dp) :: K0 = 1.7d0 ! Kolmogorov constant

      real(dp) :: bisection_a, bisection_b, bisection_c, bisection_val_1, bisection_val_2

      include 'formats'
      if (gradr > gradL) then
         ! Convection zone

         Q = chiT/chiRho ! 'Q' param  C&G 14.24
         if (MLT_option == 'Cox' .or. MLT_option == 'TDC' .or. MLT_option == 'CM1991' .or. MLT_option == 'CGM1996') then ! this assumes optically thick
            a0 = 9d0/4d0
            convective_conductivity = &
               Cp*grav*pow2(Lambda)*rho*(sqrt(Q*rho/(2d0*P)))/9d0 ! erg / (K cm sec)
            radiative_conductivity = &
               (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
            A = convective_conductivity / radiative_conductivity !  unitless.
            chi = radiative_conductivity / (Cp * rho)
         else
            select case(trim(MLT_option))
            case ('Henyey')
               ff1=1.0d0/Henyey_MLT_nu_param
               ff2=0.5d0 
               ff3=8.0d0/Henyey_MLT_y_param
               ff4=1.0d0/Henyey_MLT_y_param
            case ('ML1')
               ff1=0.125d0 
               ff2=0.5d0 
               ff3=24.0d0
               ff4=0.0d0
            case ('ML2')
               ff1=1.0d0
               ff2=2.0d0
               ff3=16.0d0
               ff4=0.0d0
            case ('Mihalas')
               ff1=0.125d0 
               ff2=0.5d0 
               ff3=16.0d0
               ff4=2.0d0
            case default
               write(*,'(3a)') 'Error: ', trim(MLT_option), ' is not an allowed MLT option'
               call mesa_error(__FILE__,__LINE__)
            end select

            omega = Lambda*rho*opacity
            ff4_omega2_plus_1 = ff4/pow2(omega) + 1d0
            a0 = (3d0/16d0)*ff2*ff3/ff4_omega2_plus_1
            A_1 = 4d0*Cp*sqrt(ff1*P*Q*rho)
            A_2 = mixing_length_alpha*omega*ff4_omega2_plus_1
            A_numerator = A_1*A_2
            A_denom = ff3*crad*clight*pow3(T)
            A = A_numerator/A_denom   
         end if 

         ! 'B' param  C&G 14.81
         Bcubed = (pow2(A)/a0)*(gradr - gradL)   

         ! now solve cubic equation for convective efficiency, Gamma
         ! a0*Gamma^3 + Gamma^2 + Gamma - a0*Bcubed == 0   C&G 14.82,
         ! leave it to Mathematica to find an expression for the root we want      
         delta = a0*Bcubed               
         f = -2d0 + 9d0*a0 + 27d0*a0*a0*delta
         if (f > 1d100) then
            f0 = f
         else
            f0 = pow2(f) + 4d0*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)
            if (f0 <= 0d0) then
               f0 = f
            else
               f0 = sqrt(f0)         
            end if
         end if
         f1 = -2d0 + 9d0*a0 + 27d0*a0*a0*delta + f0  
         if (f1 <= 0d0) return
         f1 = pow(f1,one_third)     
         f2 = 2d0*two_13*(1d0 - 3d0*a0) / f1       
         Gamma = (four_13*f1 + f2 - 2d0) / (6d0*a0)

         ! CM and CGM Implementation
         if (MLT_option == 'CM1991' .or. MLT_option == 'CGM1996') then
            if (MLT_option == 'CM1991') then
               bisection_a = 0
               bisection_b = (sqrt(1 + 4 * delta%val) - 1) / 2
               
               ! I'm not sure if this is necessary but for small values of delta
               ! the upper bound is extremely close to the root. If due to round
               ! off error the bounds don't bracket the root then bisection cannot
               ! happen. So we just artificially increase the upper bound when we 
               ! suspect round off error can occur.
               if (bisection_b < 1) then
                  bisection_b = bisection_b + 1d-2 
               end if

               fn_val = cm1991(bisection_a) * cm1991(bisection_b)

               if (fn_val > 0) then
                  ierr = 1
                  write(*,*) 'failed to get Gamma with bisection (Bad bracket)', ierr
               end if

               do while (bisection_b - bisection_a > 1d-12)
                  bisection_c = (bisection_a + bisection_b) / 2d0

                  if (bisection_c == bisection_a .or. bisection_c == bisection_b) then
                     ! For some reason the bisection is unable to improve in certain places.
                     ! It is most likely a floating point precision problem but that shouldn't
                     ! happen at real*8 unless the root is at a large enough value that the 
                     ! precision is low.
                     exit
                  end if
                  bisection_val_1 = cm1991(bisection_a)
                  bisection_val_2 = cm1991(bisection_c)
                  if (bisection_val_1 == 0) then
                     bisection_b = bisection_b
                  else if (bisection_val_2 == 0) then
                     bisection_a = bisection_b
                  else if (bisection_val_1 * bisection_val_2 < 0) then
                     bisection_b = bisection_c
                  else 
                     bisection_a = bisection_c
                  end if
                  gamma_val = (bisection_a + bisection_b) / 2d0
               end do

               fn_val = cm1991(gamma_val)
               ierr = 0
               if (abs(fn_val) > 1d-6) then 
                  ierr = 2
               end if
            else if (MLT_option == 'CGM1996') then
               bisection_a = 0
               bisection_b = (sqrt(1 + 4 * delta%val) - 1) / 2

               if (bisection_b < 1) then
                  bisection_b = bisection_b + 1d-2 
               end if
               
               fn_val = cgm1996(bisection_a) * cgm1996(bisection_b)

               if (fn_val > 0) then
                  ierr = 1
                  write(*,*) 'failed to get Gamma with bisection (Bad bracket)', ierr
               end if

               do while (bisection_b - bisection_a > 1d-12)
                  bisection_c = (bisection_a + bisection_b) / 2d0

                  if (bisection_c == bisection_a .or. bisection_c == bisection_b) then
                     ! For some reason the bisection is unable to improve in certain places.
                     ! It is most likely a floating point precision problem but that shouldn't
                     ! happen at real*8 unless the root is at a large enough value that the 
                     ! precision is low.
                     exit
                  end if
                  bisection_val_1 = cgm1996(bisection_a)
                  bisection_val_2 = cgm1996(bisection_c)
                  if (bisection_val_1 == 0) then
                     bisection_b = bisection_b
                  else if (bisection_val_2 == 0) then
                     bisection_a = bisection_b
                  else if (bisection_val_1 * bisection_val_2 < 0) then
                     bisection_b = bisection_c
                  else 
                     bisection_a = bisection_c
                  end if
                  gamma_val = (bisection_a + bisection_b) / 2d0
               end do

               fn_val = cgm1996(gamma_val)
               ierr = 0
               if (abs(fn_val) > 1d-6) then 
                  ierr = 2
               end if
            end if

            
            if (ierr /= 0 .or. is_bad(gamma_val)) then
               if (ierr == 2) then
                  if (MLT_option == "CM1991") then
                     fn_val = cm1991(gamma_val)
                  else
                     fn_val = cgm1996(gamma_val)
                  end if
                  ierr = 2
                  write(*,*) 'failed to get Gamma with bisection (Bad root)', ierr, gamma_val, fn_val
               end if
               return
            end if

            ! CM (8)
            Sigma = 2d0*gamma_val+1d0
            Sigma = Sigma*Sigma - 1d0
            dSigma_dGamma = 8d0*gamma_val + 4d0

            
            if (MLT_option == 'CM1991') then
               OmegaCM = Phi_CM(Sigma)/Phi_MLT(Sigma)
               dOmegaCM_dGamma = dSigma_dGamma*OmegaCM*(dPhi_CM_dSigma(Sigma)/Phi_CM(Sigma) &
                     - dPhi_MLT_dSigma(Sigma)/Phi_MLT(Sigma))

               Gamma%val = gamma_val
               Gamma%d1Array(1:33) = delta%d1Array(1:33)/(3d0*a0%val*gamma_val*gamma_val*OmegaCM + a0%val &
               *gamma_val*gamma_val*gamma_val*dOmegaCM_dGamma + 2d0*gamma_val + 1d0)
            else if (MLT_option == 'CGM1996') then
               OmegaCGM = Phi_CGM(Sigma)/Phi_MLT(Sigma)
               dOmegaCGM_dGamma = dSigma_dGamma*OmegaCGM*(dPhi_CGM_dSigma(Sigma)/Phi_CGM(Sigma) &
                     - dPhi_MLT_dSigma(Sigma)/Phi_MLT(Sigma))

               Gamma%val = gamma_val
               Gamma%d1Array(1:33) = delta%d1Array(1:33)/(3d0*a0%val*gamma_val*gamma_val*OmegaCGM + a0%val &
               *gamma_val*gamma_val*gamma_val*dOmegaCGM_dGamma + 2d0*gamma_val + 1d0)
            end if
         end if

         
         if (Gamma <= 0d0) return

         if (MLT_option == 'CM1991' .or. MLT_option == 'CGM1996') then
            Zeta = Phi_CM(Sigma)
            Zeta = Zeta/(1d0 + Zeta)

            Zeta%d1Array(:) = 1d0/(1d0 + Phi_CM(Sigma))**2*dPhi_CM_dSigma(Sigma) &
                   *dSigma_dGamma*Gamma%d1Array(:)

            ! CM (7)
            S = (81d0 / 2d0) * Sigma
            if (MLT_option == 'CGM1996') then
               ! average convection velocity   GGM 1996 (89)
               conv_vel = chi / Lambda * sqrt(cgm_FF3(S) * cgm_FF4(S))
            else if (MLT_option == 'CM1991') then
               ! average convection velocity   GM 1991 (41)
               conv_vel = chi * sqrt(2d0 * S) / Lambda
            end if 
            D = conv_vel*Lambda/3d0 ! TODO: Should this be different?
         else 
            ! average convection velocity   C&G 14.86b
            conv_vel = mixing_length_alpha*sqrt(Q*P/(8d0*rho))*Gamma / A
            D = conv_vel*Lambda/3d0     ! diffusion coefficient [cm^2/sec]

            !Zeta = pow3(Gamma)/Bcubed  ! C&G 14.80
            Zeta = exp(3d0*log(Gamma) - log(Bcubed)) ! write it this way to avoid overflow problems
         end if
      else
         ! Radiative zone, because this means that gradr < gradL
         Gamma = -1d99
         Zeta = 0d0
         conv_vel = 0d0
         D = 0d0
      end if

      ! Zeta must be >= 0 and <= 1.
      ! By construction (above) it cannot be less than zero,
      ! so we just check that it is a valid number and is not greater
      ! than one.
      if (is_bad(Zeta%val)) return
      if (Zeta > 1d0) then
         Zeta = 1d0
      end if

      gradT = (1d0 - Zeta)*gradr + Zeta*gradL  ! C&G 14.79
      Y_face = gradT - gradL

      if (Y_face > 0d0) then
         mixing_type = convective_mixing
      end if


      contains

      real(dp) function cgm_FF3(cgmS)
         real(dp), intent(in) :: cgmS
         real(dp), parameter:: cgm_a = 0.00101392d0, cgm_b = 0.000017848d0
         
         cgm_FF3 = pow3(K0 / 1.5) * cgm_a * pow2(cgmS) / (1 + sqrt(1 + cgm_b * pow2(cgmS)))
      end function cgm_FF3


      real(dp) function cgm_FF4(cgmS)
         real(dp), intent(in) :: cgmS
         real(dp), parameter:: cgm_c = 6.39899d0, cgm_d = 2.256815d0, &
                               cgm_e = 0.000777055d0, cgm_f = 0.868589d0
         
         cgm_FF4 = cgm_c + (cgm_d * (-1 + cgm_e *  pow(cgmS, cgm_f))) / (1 + cgm_e * pow(cgmS, cgm_f))
      end function cgm_FF4

      real(dp) function cm1991(xx)
         use const_def, only: dp
         real(dp), intent(in) :: xx

         real(dp) :: cm_Sigma, cm_OmegaCM

         ! write(*,*) 'calling CM1991...'

         ! This is only required for the first step of bisection.
         if (xx == 0) then
            cm1991 = -1
            return
         end if

         cm_Sigma = (2d0*xx+1d0)**2-1d0
         cm_OmegaCM = Phi_CM(cm_Sigma)/Phi_MLT(cm_Sigma)

         if (abs(delta%val) > 1) then
            cm1991 = (cm_OmegaCM*a0%val*xx**3 + xx**2 + xx) / delta%val - 1
         else 
            cm1991 = (cm_OmegaCM*a0%val*xx**3 + xx**2 + xx) - delta%val
         end if
      end function cm1991

      
      real(dp) function Phi_CM(cm_Sigma)
         use const_def, only: dp
         real(dp), intent(in) :: cm_Sigma
         real(dp), parameter :: aa1 = 24.868d0 ! constants from CM91, eqn (33)
         real(dp), parameter :: aa2 = 0.097666d0
         real(dp), parameter :: emm = 0.14972d0
         real(dp), parameter :: enn = 0.18931d0
         real(dp), parameter :: pee = 1.8503d0

         Phi_CM = aa1*cm_Sigma**emm*pow(pow(1d0+aa2*cm_Sigma, enn)-1d0, pee) ! CM91, eqn (32)

      end function Phi_CM

      real(dp) function dPhi_CM_dSigma(cm_Sigma)
         use const_def, only: dp
         real(dp), intent(in) :: cm_Sigma
         real(dp), parameter :: aa1 = 24.868d0 ! constants from CM91, eqn (33)
         real(dp), parameter :: aa2 = 0.097666d0
         real(dp), parameter :: emm = 0.14972d0
         real(dp), parameter :: enn = 0.18931d0
         real(dp), parameter :: pee = 1.8503d0
         ! real(dp), parameter :: share = 0.1d0

         ! dPhi_CM_dSigma = Phi_CM(Sigma)*(emm/Sigma + pee*enn*aa2*(1d0+aa2*Sigma)**(enn-1d0) &
         !      /((1d0+aa2*Sigma)**enn-1d0))
         dPhi_CM_dSigma = Phi_CM(cm_Sigma)*(emm/cm_Sigma + pee*enn*aa2*pow(1d0+aa2*cm_Sigma, enn-1d0) &
               /(pow(1d0+aa2*cm_Sigma, enn)-1d0))
      end function dPhi_CM_dSigma

      real(dp) function Phi_MLT(Sigma_MLT)
         use const_def, only: dp
         real(dp), intent(in) :: Sigma_MLT

         Phi_MLT = 0.5*a0%val/Sigma_MLT*(sqrt(1d0+Sigma_MLT)-1d0)**3
         ! Phi_MLT = 0.5*a0/Sigma*pow(sqrt(1d0+Sigma)-1d0, 3)
      end function Phi_MLT

      real(dp) function dPhi_MLT_dSigma(Sigma_MLT)
         use const_def, only: dp
         real(dp), intent(in) :: Sigma_MLT

         dPhi_MLT_dSigma = Phi_MLT(Sigma_MLT)*(1.5d0/(1d0+Sigma_MLT-sqrt(1d0+Sigma_MLT)) - 1d0/Sigma_MLT)
      end function dPhi_MLT_dSigma

      
      real(dp) function cgm1996(xx)
         ! returns with ierr = 0 if was able to evaluate f and df/dx at x
         ! if df/dx not available, it is okay to set it to 0
         use const_def, only: dp
         real(dp), intent(in) :: xx

         ! real(dp) :: Sigma, Phi_CM, Phi_MLT, OmegaCM
         real(dp) :: cgm_Sigma, cgm_OmegaCGM


         ! This is only required for the first step of bisection.
         if (xx == 0) then
            cgm1996 = -1
            return
         end if

         cgm_Sigma = (2d0*xx+1d0)**2-1d0
         cgm_OmegaCGM = Phi_CGM(cgm_Sigma)/Phi_MLT(cgm_Sigma)  


         if (abs(delta%val) > 1) then
            cgm1996 = (cgm_OmegaCGM*a0%val*xx**3 + xx**2 + xx) / delta%val - 1
         else 
            cgm1996 = (cgm_OmegaCGM*a0%val*xx**3 + xx**2 + xx) - delta%val
         end if

      end function cgm1996

      
         real(dp) function Phi_CGM(cgm_Sigma)
           use const_def, only: dp
           real(dp), intent(in) :: cgm_Sigma

           Phi_CGM = cgm_FF1(cgm_Sigma)*cgm_FF2(cgm_Sigma)
         end function Phi_CGM

         real(dp) function dPhi_CGM_dSigma(cgm_Sigma)
           use const_def, only: dp
           real(dp), intent(in) :: cgm_Sigma

           dPhi_CGM_dSigma = cgm_FF1(cgm_Sigma)*dFF2_dSigma(cgm_Sigma) + cgm_FF2(cgm_Sigma)*dFF1_dSigma(cgm_Sigma)
         end function dPhi_CGM_dSigma

         real(dp) function cgm_FF1(cgm_Sigma)
           use const_def, only: dp
           real(dp), intent(in) :: cgm_Sigma
           real(dp), parameter :: aa = 10.8654d0
           real(dp), parameter :: bb = 0.00489073d0
           real(dp), parameter :: kay = 0.149888d0
           real(dp), parameter :: emm = 0.189238d0
           real(dp), parameter :: enn = 1.85011

           cgm_FF1 = pow(K0/1.5, 3d0)*aa*pow(cgm_Sigma, kay)*pow(pow(1d0+bb*cgm_Sigma, emm)-1d0, enn)
         end function cgm_FF1

         real(dp) function dFF1_dSigma(cgm_Sigma)
           use const_def, only: dp
           real(dp), intent(in) :: cgm_Sigma
           real(dp), parameter :: aa = 10.8654d0
           real(dp), parameter :: bb = 0.00489073d0
           real(dp), parameter :: kay = 0.149888d0
           real(dp), parameter :: emm = 0.189238d0
           real(dp), parameter :: enn = 1.85011

           dFF1_dSigma = cgm_FF1(cgm_Sigma)*(kay/cgm_Sigma + &
                bb*emm*enn*pow(1d0+bb*cgm_Sigma, emm-1)/(pow(1d0+bb*cgm_Sigma, emm)-1d0))
         end function dFF1_dSigma

         real(dp) function cgm_FF2(cgm_Sigma)
           use const_def, only: dp
           real(dp), intent(in) :: cgm_Sigma
           real(dp), parameter :: cc = 0.0108071d0
           real(dp), parameter :: dd = 0.00301208d0
           real(dp), parameter :: ee = 0.000334441d0
           real(dp), parameter :: ff = 0.000125d0
           real(dp), parameter :: pee = 0.72d0
           real(dp), parameter :: que = 0.92d0
           real(dp), parameter :: are = 1.2d0
           real(dp), parameter :: tee = 1.5d0

           cgm_FF2 = 1d0 + cc*pow(cgm_Sigma, pee)/(1d0+dd*pow(cgm_Sigma, que)) &
                + ee*pow(cgm_Sigma, are)/(1d0+ff*pow(cgm_Sigma, tee))
         end function cgm_FF2

         real(dp) function dFF2_dSigma(cgm_Sigma)
           use const_def, only: dp
           real(dp), intent(in) :: cgm_Sigma
           real(dp), parameter :: cc = 0.0108071d0
           real(dp), parameter :: dd = 0.00301208d0
           real(dp), parameter :: ee = 0.000334441d0
           real(dp), parameter :: ff = 0.000125d0
           real(dp), parameter :: pee = 0.72d0
           real(dp), parameter :: que = 0.92d0
           real(dp), parameter :: are = 1.2d0
           real(dp), parameter :: tee = 1.5d0

           dFF2_dSigma = cc*pow(cgm_Sigma, pee-1d0)*(pee*(1d0+dd*pow(cgm_Sigma, que))-dd*que*pow(cgm_Sigma, que)) &
                /(1d0+dd*pow(cgm_Sigma, que))**2 &
                + ee*pow(cgm_Sigma, are-1d0)*(are*(1d0+ff*pow(cgm_Sigma, tee))-ff*tee*pow(cgm_Sigma, tee)) &
                /(1d0+ff*pow(cgm_Sigma, tee))**2
         end function dFF2_dSigma

   end subroutine calc_MLT

end module MLT
