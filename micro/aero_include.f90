!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
!
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later version.
!
! This software is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this
! code; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================

Module aero_include

implicit none

!salt source nudging time scale
real, parameter :: time_salt= 3600.

!variables for Ginoux et al.(2001) dust sources
integer, parameter :: nx_source=360,ny_source=180
real, dimension(nx_source,ny_source) :: source
real, dimension(nx_source) :: lon_source
real, dimension(ny_source) :: lat_source

!----------------------------Roughness data for deposition-----------------------
! number of land use class and season
integer, parameter :: nc = 15,ns=5 

! Z0 roughness length, A: radius of collectors, alpha, gamma
! other constants
real, dimension(nc, ns):: Z0, A, alpha, ggamma  
              
! viscosity of water at different temperature 
integer, parameter :: nt = 12 

! lookup table size for saving droplet falling speed.
integer, parameter:: nsize = 34
       
! temperature, tt, water dynamic viscosity dvcosity, N.s.m-2	 
real, dimension(nt):: tt, dvcosity
      
! data ecosystem 
       
data z0(1:nc,1) /0.80,2.65,0.85,1.05,1.15,0.1,0.1,0.04,&
                0.03,0.10,0.03,0.01,0.00,0.00,1.0/
data z0(1:nc,2) /0.80,2.65,0.85,1.05,1.15,0.1,0.1,0.04, & 
                0.03,0.10,0.03,0.01,0.00,0.00,1.0/
data z0(1:nc,3) /0.80,2.65,0.85,1.05,1.15,0.1,0.1,0.04, & 
                0.03,0.10,0.03,0.01,0.00,0.00,1.0/
data z0(1:nc,4) /0.80,2.65,0.85,1.05,1.15,0.1,0.1,0.04, & 
                0.03,0.10,0.03,0.01,0.00,0.00,1.0/
data z0(1:nc,5) /0.80,2.65,0.85,1.05,1.15,0.1,0.1,0.04, & 
                0.03,0.10,0.03,0.01,0.00,0.00,1.0/
       
data A(1:nc, 1) /2.0,5.0,2.0, 5.0,5.0,2.0,2.0,0.0,0.0, &
               10.0,10.0,0.0,0.0,0.0,10.0/
data A(1:nc, 2) /2.0,5.0,2.0, 5.0,5.0,2.0,2.0,0.0,0.0, &
               10.0,10.0,0.0,0.0,0.0,10.0/
data A(1:nc, 3) /2.0,5.0,5.0,10.0,5.0,5.0,5.0,0.0,0.0,  &
               10.0,10.0,0.0,0.0,0.0,10.0/
data A(1:nc, 4) /2.0,5.0,5.0,10.0,5.0,5.0,5.0,0.0,0.0,  &
               10.0,10.0,0.0,0.0,0.0,10.0/
data A(1:nc, 5) /2.0,5.0,2.0, 5.0,5.0,2.0,2.0,0.0,0.0,  &
               10.0,10.0,0.0,0.0,0.0,10.0/
       
! doesn't change with season       
data alpha(1:nc, 1) /1.0,  0.6,  1.1,  0.8,   0.8,   1.2, 1.2, 50.0,  &
               50.0,  1.3,  2.0, 50.0, 100.0, 100.0, 1.5/
data alpha(1:nc, 2) /1.0,  0.6,  1.1,  0.8,   0.8,   1.2, 1.2, 50.0,  &
               50.0,  1.3,  2.0, 50.0, 100.0, 100.0, 1.5/
data alpha(1:nc, 3) /1.0,  0.6,  1.1,  0.8,   0.8,   1.2, 1.2, 50.0,  &
               50.0,  1.3,  2.0, 50.0, 100.0, 100.0, 1.5/
data alpha(1:nc, 4) /1.0,  0.6,  1.1,  0.8,   0.8,   1.2, 1.2, 50.0 , &
               50.0,  1.3,  2.0, 50.0, 100.0, 100.0, 1.5/
data alpha(1:nc, 5) /1.0,  0.6,  1.1,  0.8,   0.8,   1.2, 1.2, 50.0 , &
               50.0,  1.3,  2.0, 50.0, 100.0, 100.0, 1.5/
       
! ggamma does'n change with season        
data ggamma(1:nc, 1) /0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, &
                     0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.56/ 
data ggamma(1:nc, 2) /0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, &
                     0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.56/ 
data ggamma(1:nc, 3) /0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, &
                     0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.56/ 
data ggamma(1:nc, 4) /0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, &
                     0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.56/ 
data ggamma(1:nc, 5) /0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, &
                     0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.56/ 

!Filling the data for water viscosity
!Temperature (C)
data tt(1:nt) /0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100/
!Dynamic viscosity of water (kg/m/s) x 1.e-6 (ie. 1787.e-6)
data dvcosity(1:nt) /1787,1519,1307,1002,798,653,547,467,404,355,315,282/
       
END MODULE aero_include
