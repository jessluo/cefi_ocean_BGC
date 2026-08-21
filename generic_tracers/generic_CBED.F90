module generic_CBED

   use g_tracer_utils, only : g_tracer_type, g_tracer_get_common, g_tracer_get_domain
   use g_tracer_utils, only : g_tracer_set_values, g_tracer_get_values
   use g_tracer_utils, only : g_tracer_get_pointer
   use g_tracer_utils, only : register_diag_field=>g_register_diag_field, g_send_data
   use cobalt_types,   only : generic_COBALT_type, phytoplankton, missing_value1, sperd, spery, epsln
   use cobalt_types,   only : SMALL, MEDIUM, LARGE, DIAZO, NUM_PHYTO
   use time_manager_mod,  only: time_type
   use field_manager_mod, only: fm_string_len
   use mpp_domains_mod,  only : domain2D,mpp_define_io_domain
   use data_override_mod, only: data_override
   use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, read_restart, write_restart
   use fms2_io_mod, only: register_restart_field, register_axis, register_field
   use fms_mod, only: error_mesg, NOTE, WARNING, FATAL
   use mpp_mod,           only: stdout
   use fms_mod,           only: stdout
   use FMS_co2calc_mod, only : FMS_co2calc_point
   !use, intrinsic :: ieee_arithmetic ! for checking presence of NaN or inf

   implicit none; private

   character(len=fm_string_len), parameter :: mod_name       = 'generic_CBED'
   character(len=fm_string_len), parameter :: package_name   = 'generic_cbed'

   public generic_CBED_update_from_source
   public generic_CBED_init, generic_CBED_end
   public generic_CBED_reg_diagnostics, generic_CBED_send_diagnostics

   integer, parameter :: nk_cbed = 20    ! Number of benthic layers

   type generic_CBED_type
      ! TODO: change read_porosity_from_file into a namelist variable
      logical :: read_porosity_from_file = .true.   ! flag to read porosity from file
      logical :: use_depth_dependent_OM_frac = .true.   ! flag to calculate fractions of total organic matter flux assigned to each reactivity class as a func of bathymetric depth (m)
      logical :: do_adaptive_time_stepping = .true.   ! flag to use adaptive time stepping | sub cycle dt over n steps to
      ! ensure that the change in tracer concentration in each step does not exceed a certain threshold.
      ! This is to prevent negative when reaction rates are high and the time step is too large.

      ! State variables

      real, dimension(:,:,:), allocatable :: f_o2   ! tracer o2 concentration field
      real, dimension(:,:,:), allocatable :: f_om1   ! tracer organic matter 1 (fast reacting) concentration field
      real, dimension(:,:,:), allocatable :: f_om2   ! tracer organic matter 2 (medium reacting) concentration field
      real, dimension(:,:,:), allocatable :: f_om3   ! tracer organic matter 3 (slow reacting) concentration field
      real, dimension(:,:,:), allocatable :: f_nh4   ! tracer nh4 (ammonium) concentration field
      real, dimension(:,:,:), allocatable :: f_no3   ! tracer no3 (nitrate) concentration field
      real, dimension(:,:,:), allocatable :: f_dic   ! tracer dic (dissolved inorganic carbon) concentration field
      real, dimension(:,:,:), allocatable :: f_odu   ! tracer odu (oxygen deficit unit) concentration field
      real, dimension(:,:,:), allocatable :: f_talk   ! tracer talk (total alkalinity) concentration field
      real, dimension(:,:,:), allocatable :: f_calc   ! tracer calc (calcite) concentration field
      real, dimension(:,:,:), allocatable :: f_arag   ! tracer arag (aragonite) concentration field
      real, dimension(:,:,:), allocatable :: f_ca2   ! tracer ca2 (dissolved Ca2+ ion) concentration field
      real, dimension(:,:,:), allocatable :: f_po4   ! tracer po4 (phosphate) concentration field
      ! Diagnostics
      ! 3D diags
      real, dimension(:,:,:), allocatable :: TOC              ! total organic carbon (wt %) in sediment
      real, dimension(:,:,:), allocatable :: TIC              ! total inorganic carbon (wt %) in sediment
      real, dimension(:,:,:), allocatable :: R_om_o2          ! OM resep via aerobic process
      real, dimension(:,:,:), allocatable :: R_om_no3         ! OM resep via denitrification
      real, dimension(:,:,:), allocatable :: R_om_anaerobic   ! OM resep via other anaerobic process
      real, dimension(:,:,:), allocatable :: R_dic            ! total remineralization
      real, dimension(:,:,:), allocatable :: R_nox            ! nitrification rate
      real, dimension(:,:,:), allocatable :: R_anammox        ! anammox rate
      real, dimension(:,:,:), allocatable :: R_oduox          ! ODU reoxidation rate
      real, dimension(:,:,:), allocatable :: R_talk_org       ! TA organic
      real, dimension(:,:,:), allocatable :: R_talk_inorg     ! TA inorganic
      real, dimension(:,:,:), allocatable :: R_talk           ! total TA = TA org + TA inorg
      real, dimension(:,:,:), allocatable :: cbed_bioirri
      real, dimension(:,:,:), allocatable :: cbed_omega_calc
      real, dimension(:,:,:), allocatable :: cbed_omega_arag
      real, dimension(:,:,:), allocatable :: cbed_ph

      ! 2D diags
      real, dimension(:,:), allocatable :: o2_flux !benthic o2 flux
      real, dimension(:,:), allocatable :: nh4_flux !benthic nh4 flux
      real, dimension(:,:), allocatable :: no3_flux !benthic no3 flux
      real, dimension(:,:), allocatable :: dic_flux !benthic dic flux
      real, dimension(:,:), allocatable :: talk_flux !benthic talk flux
      real, dimension(:,:), allocatable :: odu_flux !benthic odu flux
      real, dimension(:,:), allocatable :: burial_om !organic matter burial at the bottom of sediment column
      real, dimension(:,:), allocatable :: denit
      real, dimension(:,:), allocatable :: cbed_k1
      real, dimension(:,:), allocatable :: cbed_k2
      real, dimension(:,:), allocatable :: cbed_k3
      real, dimension(:,:), allocatable :: cbed_w
      !real, dimension(:,:), allocatable :: cbed_anammox
      !real, dimension(:,:), allocatable :: cbed_o2resp
      !real, dimension(:,:), allocatable :: cbed_no3resp

      ! accumulated b_terms for passing to next timestep | used due to adaptive timestepping
      real, dimension(:,:), allocatable :: cbed_b_o2_acc
      real, dimension(:,:), allocatable :: cbed_b_dic_acc
      real, dimension(:,:), allocatable :: cbed_b_nh4_acc
      real, dimension(:,:), allocatable :: cbed_b_no3_acc
      real, dimension(:,:), allocatable :: cbed_b_alk_acc
      real, dimension(:,:), allocatable :: cbed_b_odu_acc
      real, dimension(:,:), allocatable :: cbed_b_ca2_acc
      real, dimension(:,:), allocatable :: cbed_b_po4_acc

      ! 3D diags, CBED grid
      real, dimension(:,:,:), allocatable :: dz_cbed             ! cbed grid thickness
      real, dimension(:,:,:), allocatable :: z_cbed_mid          ! cbed layer mid points

      ! 3D diags (interfaces)(nk_cbed+1)
      real, dimension(:,:,:), allocatable :: cbed_Db
      real, dimension(:,:,:), allocatable :: cbed_D_o2
      real, dimension(:,:,:), allocatable :: cbed_D_dic
      real, dimension(:,:,:), allocatable :: cbed_D_nh4
      real, dimension(:,:,:), allocatable :: cbed_D_no3
      real, dimension(:,:,:), allocatable :: cbed_D_odu
      real, dimension(:,:,:), allocatable :: cbed_por
      real, dimension(:,:,:), allocatable :: cbed_svf



      integer :: id_o2                               ! tracer o2 diagnostics id
      integer :: id_om1                              ! tracer om1 diagnostics id
      integer :: id_om2                              ! tracer om2 diagnostics id
      integer :: id_om3                              ! tracer om3 diagnostics id
      integer :: id_nh4                              ! tracer nh4 diagnostics id
      integer :: id_no3                              ! tracer no3 diagnostics id
      integer :: id_dic                              ! tracer dic diagnostics id
      integer :: id_odu                              ! tracer odu diagnostics id
      integer :: id_talk                             ! tracer talk diagnostics id
      integer :: id_calc
      integer :: id_arag
      integer :: id_ca2
      integer :: id_po4
      ! 3D diags
      integer :: id_TOC
      integer :: id_TIC
      integer :: id_R_om_o2
      integer :: id_R_om_no3
      integer :: id_R_om_anaerobic
      integer :: id_R_dic
      integer :: id_R_nox
      integer :: id_R_anammox
      integer :: id_R_oduox
      integer :: id_R_talk_org
      integer :: id_R_talk_inorg
      integer :: id_R_talk
      integer :: id_cbed_bioirri
      integer :: id_cbed_omega_calc
      integer :: id_cbed_omega_arag
      integer :: id_cbed_ph

      ! 2D diags
      integer :: id_o2_flux
      integer :: id_nh4_flux
      integer :: id_no3_flux
      integer :: id_dic_flux
      integer :: id_talk_flux
      integer :: id_odu_flux
      integer :: id_burial_om
      integer :: id_denit
      integer :: id_cbed_k1
      integer :: id_cbed_k2
      integer :: id_cbed_k3
      integer :: id_cbed_w

      integer :: id_cbed_b_o2_acc
      integer :: id_cbed_b_dic_acc
      integer :: id_cbed_b_nh4_acc
      integer :: id_cbed_b_no3_acc
      integer :: id_cbed_b_alk_acc
      integer :: id_cbed_b_odu_acc
      integer :: id_cbed_b_ca2_acc
      integer :: id_cbed_b_po4_acc

      !integer :: id_cbed_anammox
      !integer :: id_cbed_o2resp
      !integer :: id_cbed_no3resp
      ! 3D diag, CBED grid
      integer :: id_dz_cbed
      integer :: id_z_cbed_mid

      ! 3D diags (interfaces)(nk_cbed+1)
      integer :: id_cbed_Db
      integer :: id_cbed_D_o2
      integer :: id_cbed_D_dic
      integer :: id_cbed_D_nh4
      integer :: id_cbed_D_no3
      integer :: id_cbed_D_odu
      integer :: id_cbed_por
      integer :: id_cbed_svf


   end type generic_CBED_type

   type(generic_CBED_type) :: cbed


   real, parameter :: pi = acos(-1.0)


   ! grid
   ! local parameters
   real, parameter :: l_cbed = 0.20           ! length of sediment domain | sediment depth (m, 20 cm)
   real, parameter :: dz1_cbed = 0.001        ! thickness of the first layer (m). For increasing thickness
   real, parameter :: rho_s = 2.5             ! solid density (g/cm³)
   real, parameter :: Db_l = 0.08             ! bioturbation length scale (m) 8 cm.
   real, parameter :: bioirri_l = 0.018       ! bioirrigation length scale (m) 1.8 cm.

   ! sediment grid and state variables (to be allocated)
   !real, allocatable :: dz_cbed(:)                 ! sediment layer thickness (m)
   !real, allocatable :: z_cbed(:)              ! sediment depth points (m)
   real :: dz_cbed(nk_cbed)              ! thickness of each cbed layers (m)
   real :: z_cbed_int(nk_cbed+1)         ! layer interfaces (m)
   real :: z_cbed_mid(nk_cbed)           ! layer mid points (m)

   ! grid param end.

   real, dimension(:,:,:), allocatable :: por       !porosity
   real, dimension(:,:,:), allocatable :: svf       !solid volume fraction (1-porosity)

   real, dimension(:,:,:), allocatable :: w       !sedimentation rate
   real, dimension(:,:), allocatable :: Db_0    !max bioturbation rate
   real, dimension(:,:,:), allocatable :: Db    !bioturbation
   real, dimension(:,:), allocatable :: bioirri_0    !max bioirrigation rate
   real, dimension(:,:,:), allocatable :: bioirri    !bioirrigation

   real, dimension(:,:,:), allocatable :: D_o2    !diffusion coefficient for o2
   real, dimension(:,:,:), allocatable :: D_dic    !diffusion coefficient for DIC
   real, dimension(:,:,:), allocatable :: D_nh4    !diffusion coefficient for NH4
   real, dimension(:,:,:), allocatable :: D_no3    !diffusion coefficient for NO3
   real, dimension(:,:,:), allocatable :: D_odu    !diffusion coefficient for ODU (H2S)
   real, dimension(:,:,:), allocatable :: D_ca2    !diffusion coefficient for Ca2+
   real, dimension(:,:,:), allocatable :: D_po4    !diffusion coefficient for PO4

   real, dimension(:,:), allocatable :: k1       ! k1 is the rate constant for the first order reaction of OM1 decomposition
   real, dimension(:,:), allocatable :: k2       ! k2 is the rate constant for the first order reaction of OM2 decomposition
   real, dimension(:,:), allocatable :: k3       ! k3 is the rate constant for the first order reaction of OM3 decomposition




!     call grid_cbed(nk_cbed, dz_cbed, z_cbed_mid, z_cbed_int)

!    ! define uniform sediment grid
!    dz_cbed = l_cbed / real(nk_cbed)
!    z_cbed_int(1) = 0.0   !this is likely the interface. dimention of z_cbed is nk_cbed+1. z_int_cbed. might need z_mid_cbed
!    do k = 1, nk_cbed
!        z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
!    end do
!
!    z_cbed_mid(1) = dz_cbed(1)/2   ! first layer mid point
!    do k = 1, nk_cbed-1
!        z_cbed_mid(k+1) = z_cbed_mid(k) + dz_cbed(k)
!    end do


contains

!! To make increasing thickness CBED grid
! Function to find the common ratio r using bisection method
   function find_r(nk_cbed, dz_first, total_height) result(r)
      integer, intent(in) :: nk_cbed
      real, intent(in) :: dz_first, total_height
      real :: r
      real :: r_low, r_high, r_mid
      real :: power, sum_geom
      integer :: i, j

      r_low = 1.0001
      r_high = 10.0
      do i = 1, 100
         r_mid = (r_low + r_high) / 2.0
         ! Compute r_mid**nk_cbed using iterative multiplication to avoid overflow
         power = 1.0
         do j = 1, nk_cbed
            power = power * r_mid
         end do
         sum_geom = dz_first * (power - 1.0) / (r_mid - 1.0)
         if (sum_geom < total_height) then
            r_low = r_mid
         else
            r_high = r_mid
         end if
      end do
      r = r_mid
   end function find_r


! Function to calculate seawater calcium (Ca2+) conc. from salinity. mol/m3
   function fn_sw_ca2(salinity, sw_density) result(dissolved_calcium)
      real, intent(in) :: salinity, sw_density
      real :: dissolved_calcium
      real :: ca_at_35, salinity_ref, ca_ratio

      ! Standard value: ~10.28 mmol/kg at Salinity 35
      ca_at_35 = 10.28 * 1e-3 * sw_density
      salinity_ref = 35.0
      ! Calculate the ratio
      ca_ratio = ca_at_35 / salinity_ref
      ! Calculate calcium concentration for the given salinity
      dissolved_calcium = salinity * ca_ratio  !unit: mol/m3
   end function fn_sw_ca2


   subroutine generic_CBED_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk)
      integer,     intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk
      !Locals
      type(domain2D), pointer :: domain
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file

      integer :: i,j,k !for grid.
      real    :: r ! for grid

      !real,dimension(isc:iec,jsc:jec,nk_cbed)    :: cbed_tmask
      ! Make a cbed mask. Note: it seems grid_tmask(:,:,k) does not depend on k
      ! Note that grid_tmask is already on isc:iec, jsc:jec
      !do j = jsc, jec; do i = isc, iec; do k=1,nk_cbed ;
      !   cbed_tmask(i,j,k) = grid_tmask(i,j,nk) ; enddo; enddo; enddo

      !Allocate and initialize CBED arrays for tracer concentrations and other workarrays
      allocate(cbed%f_o2(isd:ied,jsd:jed,nk_cbed));cbed%f_o2=0.0
      allocate(cbed%f_om1(isd:ied,jsd:jed,nk_cbed));cbed%f_om1=0.0
      allocate(cbed%f_om2(isd:ied,jsd:jed,nk_cbed));cbed%f_om2=0.0
      allocate(cbed%f_om3(isd:ied,jsd:jed,nk_cbed));cbed%f_om3=0.0
      allocate(cbed%f_nh4(isd:ied,jsd:jed,nk_cbed));cbed%f_nh4=0.0
      allocate(cbed%f_no3(isd:ied,jsd:jed,nk_cbed));cbed%f_no3=0.0
      allocate(cbed%f_dic(isd:ied,jsd:jed,nk_cbed));cbed%f_dic=0.0
      allocate(cbed%f_odu(isd:ied,jsd:jed,nk_cbed));cbed%f_odu=0.0
      allocate(cbed%f_talk(isd:ied,jsd:jed,nk_cbed));cbed%f_talk=0.0
      allocate(cbed%f_calc(isd:ied,jsd:jed,nk_cbed));cbed%f_calc=0.0
      allocate(cbed%f_arag(isd:ied,jsd:jed,nk_cbed));cbed%f_arag=0.0
      allocate(cbed%f_ca2(isd:ied,jsd:jed,nk_cbed));cbed%f_ca2=0.0
      allocate(cbed%f_po4(isd:ied,jsd:jed,nk_cbed));cbed%f_po4=0.0
      !Diagnostics
      ! 3D diags
      allocate(cbed%TOC(isd:ied,jsd:jed,nk_cbed));cbed%TOC=0.0
      allocate(cbed%TIC(isd:ied,jsd:jed,nk_cbed));cbed%TIC=0.0
      allocate(cbed%R_om_o2(isd:ied,jsd:jed,nk_cbed));cbed%R_om_o2=0.0
      allocate(cbed%R_om_no3(isd:ied,jsd:jed,nk_cbed));cbed%R_om_no3=0.0
      allocate(cbed%R_om_anaerobic(isd:ied,jsd:jed,nk_cbed));cbed%R_om_anaerobic=0.0
      allocate(cbed%R_dic(isd:ied,jsd:jed,nk_cbed));cbed%R_dic=0.0
      allocate(cbed%R_nox(isd:ied,jsd:jed,nk_cbed));cbed%R_nox=0.0
      allocate(cbed%R_anammox(isd:ied,jsd:jed,nk_cbed));cbed%R_anammox=0.0
      allocate(cbed%R_oduox(isd:ied,jsd:jed,nk_cbed));cbed%R_oduox=0.0
      allocate(cbed%R_talk_org(isd:ied,jsd:jed,nk_cbed));cbed%R_talk_org=0.0
      allocate(cbed%R_talk_inorg(isd:ied,jsd:jed,nk_cbed));cbed%R_talk_inorg=0.0
      allocate(cbed%R_talk(isd:ied,jsd:jed,nk_cbed));cbed%R_talk=0.0
      allocate(cbed%cbed_bioirri(isd:ied,jsd:jed,nk_cbed));cbed%cbed_bioirri=0.0
      allocate(cbed%cbed_omega_calc(isd:ied,jsd:jed,nk_cbed));cbed%cbed_omega_calc=0.0
      allocate(cbed%cbed_omega_arag(isd:ied,jsd:jed,nk_cbed));cbed%cbed_omega_arag=0.0
      allocate(cbed%cbed_ph(isd:ied,jsd:jed,nk_cbed));cbed%cbed_ph=0.0
      ! 2D diags
      allocate(cbed%o2_flux(isd:ied,jsd:jed)); cbed%o2_flux=0.0
      allocate(cbed%nh4_flux(isd:ied,jsd:jed)); cbed%nh4_flux=0.0
      allocate(cbed%no3_flux(isd:ied,jsd:jed)); cbed%no3_flux=0.0
      allocate(cbed%dic_flux(isd:ied,jsd:jed)); cbed%dic_flux=0.0
      allocate(cbed%talk_flux(isd:ied,jsd:jed)); cbed%talk_flux=0.0
      allocate(cbed%odu_flux(isd:ied,jsd:jed)); cbed%odu_flux=0.0
      allocate(cbed%burial_om(isd:ied,jsd:jed));cbed%burial_om=0.0
      allocate(cbed%denit(isd:ied,jsd:jed));cbed%denit=0.0
      allocate(cbed%cbed_k1(isd:ied,jsd:jed));cbed%cbed_k1=0.0
      allocate(cbed%cbed_k2(isd:ied,jsd:jed));cbed%cbed_k2=0.0
      allocate(cbed%cbed_k3(isd:ied,jsd:jed));cbed%cbed_k3=0.0
      allocate(cbed%cbed_w(isd:ied,jsd:jed));cbed%cbed_w=0.0

      allocate(cbed%cbed_b_o2_acc(isd:ied,jsd:jed));cbed%cbed_b_o2_acc=0.0
      allocate(cbed%cbed_b_dic_acc(isd:ied,jsd:jed));cbed%cbed_b_dic_acc=0.0
      allocate(cbed%cbed_b_nh4_acc(isd:ied,jsd:jed));cbed%cbed_b_nh4_acc=0.0
      allocate(cbed%cbed_b_no3_acc(isd:ied,jsd:jed));cbed%cbed_b_no3_acc=0.0
      allocate(cbed%cbed_b_alk_acc(isd:ied,jsd:jed));cbed%cbed_b_alk_acc=0.0
      allocate(cbed%cbed_b_odu_acc(isd:ied,jsd:jed));cbed%cbed_b_odu_acc=0.0
      allocate(cbed%cbed_b_ca2_acc(isd:ied,jsd:jed));cbed%cbed_b_ca2_acc=0.0
      allocate(cbed%cbed_b_po4_acc(isd:ied,jsd:jed));cbed%cbed_b_po4_acc=0.0

      !allocate(cbed%cbed_anammox(isd:ied,jsd:jed));cbed%cbed_anammox=0.0
      !allocate(cbed%cbed_o2resp(isd:ied,jsd:jed));cbed%cbed_o2resp=0.0
      !allocate(cbed%cbed_no3resp(isd:ied,jsd:jed));cbed%cbed_no3resp=0.0
      ! 3D diag, CBED grid
      allocate(cbed%dz_cbed(isd:ied,jsd:jed,nk_cbed));cbed%dz_cbed=0.0
      allocate(cbed%z_cbed_mid(isd:ied,jsd:jed,nk_cbed));cbed%z_cbed_mid=0.0

      ! 3D diags (interfaces)(nk_cbed+1)
      allocate(cbed%cbed_Db(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_Db=0.0
      allocate(cbed%cbed_D_o2(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_o2=0.0
      allocate(cbed%cbed_D_dic(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_dic =0.0
      allocate(cbed%cbed_D_nh4(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_nh4=0.0
      allocate(cbed%cbed_D_no3(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_no3=0.0
      allocate(cbed%cbed_D_odu(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_odu=0.0
      allocate(cbed%cbed_por(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_por=0.0
      allocate(cbed%cbed_svf(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_svf=0.0


      allocate(por(isd:ied,jsd:jed,nk_cbed+1));        por=0.0  !initialize to zero
      allocate(svf(isd:ied,jsd:jed,nk_cbed+1));        svf=0.0  !solid volume fraction

      allocate(w(isd:ied,jsd:jed,nk_cbed+1));        w=0.0      !adding sedimentation rate initalize
      allocate(Db_0(isd:ied,jsd:jed));               Db_0=0.0   !bioturbation_0 init.
      allocate(Db(isd:ied,jsd:jed,nk_cbed+1));       Db=0.0     !bioturbation init.
      allocate(bioirri_0(isd:ied,jsd:jed));          bioirri_0=0.0   !bioirrigation_0 init.
      allocate(bioirri(isd:ied,jsd:jed,nk_cbed));    bioirri=0.0     !bioturbation init.

      allocate(D_o2(isd:ied,jsd:jed,nk_cbed+1)); D_o2=0.0     ! D_o2 init.
      allocate(D_dic(isd:ied,jsd:jed,nk_cbed+1)); D_dic=0.0
      allocate(D_nh4(isd:ied,jsd:jed,nk_cbed+1)); D_nh4=0.0
      allocate(D_no3(isd:ied,jsd:jed,nk_cbed+1)); D_no3=0.0
      allocate(D_odu(isd:ied,jsd:jed,nk_cbed+1)); D_odu=0.0
      allocate(D_ca2(isd:ied,jsd:jed,nk_cbed+1)); D_ca2=0.0
      allocate(D_po4(isd:ied,jsd:jed,nk_cbed+1)); D_po4=0.0

      allocate(k1(isd:ied,jsd:jed)); k1=0.0
      allocate(k2(isd:ied,jsd:jed)); k2=0.0
      allocate(k3(isd:ied,jsd:jed)); k3=0.0

      !if (cbed%read_porosity_from_file) then
      !   call data_override('OCN', 'por', por(isc:iec,jsc:jec,nk_cbed+1), model_time, override=.true.)
      !   svf(isc:iec,jsc:jec,nk_cbed+1) = 1.0 - por(isc:iec,jsc:jec,nk_cbed+1)
      !else
      !   por(isc:iec,jsc:jec,nk_cbed+1) = 0.8
      !   svf(isc:iec,jsc:jec,nk_cbed+1) = 0.2
      !endif


      ! Grid does not change with time, so can be define only once.

      !!!! define uniform sediment grid
      !dz_cbed = l_cbed / real(nk_cbed)
      !z_cbed_int(1) = 0.0   !this is likely the interface. dimention of z_cbed is nk_cbed+1. z_int_cbed. might need z_mid_cbed
      !do k = 1, nk_cbed
      !   z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
      !end do
      !
      !z_cbed_mid(1) = dz_cbed(1)/2   ! first layer mid point
      !do k = 1, nk_cbed-1
      !   z_cbed_mid(k+1) = z_cbed_mid(k) + dz_cbed(k)
      !end do


      ! generate grid automatically according to supplied values (length of sediment column, number of layers and thicness of the first layer)
      ! Compute the common ratio r
      r = find_r(nk_cbed, dz1_cbed, l_cbed)

      ! Generate dz_cbed as geometric progression
      dz_cbed(1) = dz1_cbed
      do k = 2, nk_cbed
         dz_cbed(k) = dz_cbed(k-1) * r
      end do

      ! Optional: Scale to ensure exact total sum (due to floating-point precision)
      ! real :: actual_sum
      ! actual_sum = sum(dz_cbed)
      ! if (abs(actual_sum - l_cbed) > 1e-10) then
      !   dz_cbed = dz_cbed * l_cbed / actual_sum
      ! end if

      ! Compute layer interfaces
      z_cbed_int(1) = 0.0
      do k = 1, nk_cbed
         z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
      end do

      ! Compute layer midpoints
      do k = 1, nk_cbed
         z_cbed_mid(k) = (z_cbed_int(k) + z_cbed_int(k+1)) / 2.0
      end do


   end subroutine generic_CBED_init

   subroutine generic_CBED_reg_diagnostics(axes,init_time)
      USE diag_manager_mod, ONLY: register_diag_field, diag_axis_init
      integer,         intent(in) :: axes(3)
      type(time_type), intent(in) :: init_time
      !Locals
      integer :: k, id_layer, id_layer_i
      real :: cbed_layers(1:nk_cbed)
      real :: cbed_layers_i(1:nk_cbed+1)

      !!BEGIN read_restart code block
      !Niki: This code block does not seem to belong here and should be in the _init routine instead.
      !      But the problem with that is due to MOM6 code flow, when generic_CBED_init is called
      !      the MOM "domain" is not yet created/updated and we cannot access it, which is needed for reading the restart here.
      type(domain2D), pointer :: domain
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file
      !Resgister restarts
      restart_file = 'INPUT/generic_CBED.res.nc'
      call g_tracer_get_domain(domain)
      file_open_success=open_file(fileobj, trim(restart_file),"read", domain, is_restart=.true.)
      if (file_open_success) then
         call register_axis(fileobj,'x','x')
         call register_axis(fileobj,'y','y')
         !!< Register the domain decomposed dimensions as variables so that the combiner can work correctly
         !call register_field(fileobj, "x", "double", (/"x"/))
         !call register_field(fileobj, "y", "double", (/"y"/))
         call register_axis(fileobj,'lev',nk_cbed)
         ! register the restart variables
         call register_restart_field(fileobj, "cbed_o2", cbed%f_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om1", cbed%f_om1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om2", cbed%f_om2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om3", cbed%f_om3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_nh4", cbed%f_nh4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_no3", cbed%f_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_dic", cbed%f_dic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_odu", cbed%f_odu, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_talk", cbed%f_talk, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_calc", cbed%f_calc, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_arag", cbed%f_arag, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_ca2", cbed%f_ca2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_po4", cbed%f_po4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_b_o2_acc", cbed%cbed_b_o2_acc, (/"x","y"/))     ! b_o2_acc etc terms are needed in restart bcz they are passed to next time step. If not as restart, they will loose information when extending a run (recheck reasoning)
         call register_restart_field(fileobj, "cbed_b_dic_acc", cbed%cbed_b_dic_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_nh4_acc", cbed%cbed_b_nh4_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_no3_acc", cbed%cbed_b_no3_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_alk_acc", cbed%cbed_b_alk_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_odu_acc", cbed%cbed_b_odu_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_ca2_acc", cbed%cbed_b_ca2_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_po4_acc", cbed%cbed_b_po4_acc, (/"x","y"/))

         call read_restart(fileobj)
      endif
      !!END read_restart code block

      !!Register diagnostics
      !Niki: I am unsure of the diag axis thingy, ask Yi-Cheng
      !Define cbed layer axis, the x,y axes are the same as MOM6 since the horizontal grids are the same
      do k=1,nk_cbed; cbed_layers(k) = k; enddo
      id_layer = diag_axis_init('cbedlayer', cbed_layers, 'None', 'z', long_name='Benthos Layer', direction=-1)
      do k=1,nk_cbed+1; cbed_layers_i(k) = k; enddo
      id_layer_i = diag_axis_init('cbedlayer_i', cbed_layers_i, 'None', 'z', long_name='Benthos Layer Interface', direction=-1)

      cbed%id_o2 = register_diag_field(package_name, 'cbed_o2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed oxygen concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_om1 = register_diag_field(package_name, 'cbed_om1_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM1 concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_om2 = register_diag_field(package_name, 'cbed_om2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM2 concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_om3 = register_diag_field(package_name, 'cbed_om3_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM3 concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_nh4 = register_diag_field(package_name, 'cbed_nh4_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed ammonium concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_no3 = register_diag_field(package_name, 'cbed_no3_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed nitrate concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_dic = register_diag_field(package_name, 'cbed_dic_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed DIC concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_odu = register_diag_field(package_name, 'cbed_odu_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed ODU concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_talk = register_diag_field(package_name, 'cbed_talk_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed total alkalinity concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_calc = register_diag_field(package_name, 'cbed_calc_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'calcite concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_arag = register_diag_field(package_name, 'cbed_arag_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'aragonite concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_ca2 = register_diag_field(package_name, 'cbed_ca2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed Ca2+ concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_po4 = register_diag_field(package_name, 'cbed_po4_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM1 concentration', 'mol m-3', missing_value = missing_value1)
      ! diags
      ! 3D diags
      cbed%id_TOC = register_diag_field(package_name, 'cbed_TOC', (/axes(1),axes(2),id_layer/), init_time,&
         'Total Organic Carbon in sediment', 'wt %', missing_value = missing_value1)
      cbed%id_TIC = register_diag_field(package_name, 'cbed_TIC', (/axes(1),axes(2),id_layer/), init_time,&
         'Total Inorganic Carbon in sediment', 'wt %', missing_value = missing_value1)
      cbed%id_R_om_o2 = register_diag_field(package_name, 'cbed_R_om_o2', (/axes(1),axes(2),id_layer/), init_time,&
         'aerobic respiration in sediment 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_om_no3 = register_diag_field(package_name, 'cbed_R_om_no3', (/axes(1),axes(2),id_layer/), init_time,&
         'OM respiration via denitrification in sediment 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_om_anaerobic = register_diag_field(package_name, 'cbed_R_om_anaerobic', (/axes(1),axes(2),id_layer/), init_time,&
         'OM respiration via other anaerobic processes in sediment 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_dic = register_diag_field(package_name, 'cbed_R_dic', (/axes(1),axes(2),id_layer/), init_time,&
         'DIC produced in sediment via OM remineralization 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_nox = register_diag_field(package_name, 'cbed_R_nox', (/axes(1),axes(2),id_layer/), init_time,&
         'nitrification in sediment', 'mol m-3 s-1', missing_value = missing_value1)
      cbed%id_R_anammox = register_diag_field(package_name, 'cbed_R_anammox', (/axes(1),axes(2),id_layer/), init_time,&
         'anammox in sediment', 'mol N2 m-3 s-1', missing_value = missing_value1)
      cbed%id_R_oduox = register_diag_field(package_name, 'cbed_R_oduox', (/axes(1),axes(2),id_layer/), init_time,&
         'ODU reoxidation in sediment', 'mol m-3 s-1', missing_value = missing_value1)
      cbed%id_R_talk_org = register_diag_field(package_name, 'cbed_R_talk_org', (/axes(1),axes(2),id_layer/), init_time,&
         'net organic TA prod', 'mol m-3 s-1', missing_value = missing_value1)
      cbed%id_R_talk_inorg = register_diag_field(package_name, 'cbed_R_talk_inorg', (/axes(1),axes(2),id_layer/), init_time,&
         'net inorganic TA prod', 'mol m-3 s-1', missing_value = missing_value1)
      cbed%id_R_talk = register_diag_field(package_name, 'cbed_R_talk', (/axes(1),axes(2),id_layer/), init_time,&
         'net TA prod, org+inorg', 'mol m-3 s-1', missing_value = missing_value1)
      cbed%id_cbed_bioirri = register_diag_field(package_name, 'cbed_bioirri', (/axes(1),axes(2),id_layer/), init_time,&
         'bioirrigation coefficient', 's-1', missing_value = missing_value1)
      cbed%id_cbed_omega_calc = register_diag_field(package_name, 'cbed_omega_calc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed omega calcite', 'mol/kg', missing_value = missing_value1)
      cbed%id_cbed_omega_arag = register_diag_field(package_name, 'cbed_omega_arag', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed omega aragonite', 'mol/kg', missing_value = missing_value1)
      cbed%id_cbed_ph = register_diag_field(package_name, 'cbed_ph', (/axes(1),axes(2),id_layer/), init_time,&
         'sediment pH', 'total scale', missing_value = missing_value1)

      ! 2D diags
      cbed%id_o2_flux = register_diag_field(package_name, 'cbed_o2_flux', (/axes(1),axes(2)/), init_time,&
         'benthic O2 flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_nh4_flux = register_diag_field(package_name, 'cbed_nh4_flux', (/axes(1),axes(2)/), init_time,&
         'benthic nh4 flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_no3_flux = register_diag_field(package_name, 'cbed_no3_flux', (/axes(1),axes(2)/), init_time,&
         'benthic no3 flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_dic_flux = register_diag_field(package_name, 'cbed_dic_flux', (/axes(1),axes(2)/), init_time,&
         'benthic dic flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_talk_flux = register_diag_field(package_name, 'cbed_talk_flux', (/axes(1),axes(2)/), init_time,&
         'benthic talk flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_odu_flux = register_diag_field(package_name, 'cbed_odu_flux', (/axes(1),axes(2)/), init_time,&
         'benthic odu flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_burial_om = register_diag_field(package_name, 'cbed_burial_om', (/axes(1),axes(2)/), init_time,&
         'cbed organic carbon burial', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_denit = register_diag_field(package_name, 'cbed_denit', (/axes(1),axes(2)/), init_time,&
         'cbed total denitrification (denit + anammox)', 'mol N m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_k1 = register_diag_field(package_name, 'cbed_k1', (/axes(1),axes(2)/), init_time,&
         'OM1 decay rate constant', 's-1', missing_value = missing_value1)
      cbed%id_cbed_k2 = register_diag_field(package_name, 'cbed_k2', (/axes(1),axes(2)/), init_time,&
         'OM2 decay rate constant', 's-1', missing_value = missing_value1)
      cbed%id_cbed_k3 = register_diag_field(package_name, 'cbed_k3', (/axes(1),axes(2)/), init_time,&
         'OM3 decay rate constant', 's-1', missing_value = missing_value1)
      cbed%id_cbed_w = register_diag_field(package_name, 'cbed_w', (/axes(1),axes(2)/), init_time,&
         'sedimentation rate', 'm/s', missing_value = missing_value1)

      cbed%id_cbed_b_o2_acc = register_diag_field(package_name, 'cbed_b_o2_acc', (/axes(1),axes(2)/), init_time,&
         'b_o2 to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_dic_acc = register_diag_field(package_name, 'cbed_b_dic_acc', (/axes(1),axes(2)/), init_time,&
         'b_dic to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_nh4_acc = register_diag_field(package_name, 'cbed_b_nh4_acc', (/axes(1),axes(2)/), init_time,&
         'b_nh4 to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_no3_acc = register_diag_field(package_name, 'cbed_b_no3_acc', (/axes(1),axes(2)/), init_time,&
         'b_no3 to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_alk_acc = register_diag_field(package_name, 'cbed_b_alk_acc', (/axes(1),axes(2)/), init_time,&
         'b_alk to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_odu_acc = register_diag_field(package_name, 'cbed_b_odu_acc', (/axes(1),axes(2)/), init_time,&
         'b_odu to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_ca2_acc = register_diag_field(package_name, 'cbed_b_ca2_acc', (/axes(1),axes(2)/), init_time,&
         'b_ca to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_b_po4_acc = register_diag_field(package_name, 'cbed_b_po4_acc', (/axes(1),axes(2)/), init_time,&
         'b_po4 to pass to next timestep', 'mol m-2 s-1', missing_value = missing_value1)

      !cbed%id_cbed_anammox = register_diag_field(package_name, 'cbed_anammox', (/axes(1),axes(2)/), init_time,&
      !   'cbed anammox', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_o2resp = register_diag_field(package_name, 'cbed_o2resp', (/axes(1),axes(2)/), init_time,&
      !   'cbed OM respiration by O2', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_no3resp = register_diag_field(package_name, 'cbed_no3resp', (/axes(1),axes(2)/), init_time,&
      !   'cbed OM respiration by NO3', 'mol/m2/s', missing_value = missing_value1)
      ! 3D diags, CBED grid
      cbed%id_dz_cbed = register_diag_field(package_name, 'cbed_dz_cbed', (/axes(1),axes(2),id_layer/), init_time,&
         'CBED grid layer thickess', 'm', missing_value = missing_value1)
      cbed%id_z_cbed_mid = register_diag_field(package_name, 'cbed_z_cbed_mid', (/axes(1),axes(2),id_layer/), init_time,&
         'CBED grid layer midpoints', 'm', missing_value = missing_value1)

      ! 3D diags, CBED grid interfaces
      cbed%id_cbed_Db = register_diag_field(package_name, 'cbed_Db', (/axes(1),axes(2),id_layer_i/), init_time,&
         'bioturbation coefficient', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_o2 = register_diag_field(package_name, 'cbed_D_o2', (/axes(1),axes(2),id_layer_i/), init_time,&
         'O2 molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_dic = register_diag_field(package_name, 'cbed_D_dic', (/axes(1),axes(2),id_layer_i/), init_time,&
         'DIC molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_nh4 = register_diag_field(package_name, 'cbed_D_nh4', (/axes(1),axes(2),id_layer_i/), init_time,&
         'NH4 molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_no3 = register_diag_field(package_name, 'cbed_D_no3', (/axes(1),axes(2),id_layer_i/), init_time,&
         'NO3 molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_odu = register_diag_field(package_name, 'cbed_D_odu', (/axes(1),axes(2),id_layer_i/), init_time,&
         'ODU molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_por = register_diag_field(package_name, 'cbed_por', (/axes(1),axes(2),id_layer_i/), init_time,&
         'sediment porosity', 'dimensionless', missing_value = missing_value1)
      cbed%id_cbed_svf = register_diag_field(package_name, 'cbed_svf', (/axes(1),axes(2),id_layer_i/), init_time,&
         'sediment solid volume fraction', 'dimensionless', missing_value = missing_value1)

   end subroutine generic_CBED_reg_diagnostics

   subroutine generic_CBED_send_diagnostics(model_time, isc,iec,jsc,jec, isd,ied,jsd,jed,nk, grid_tmask)
      USE diag_manager_mod, ONLY: send_data
      type(time_type),          intent(in) :: model_time
      real, dimension(:,:,:),    pointer   :: grid_tmask
      integer,                  intent(in) :: isc,iec,jsc,jec, isd,ied,jsd,jed,nk
      ! local
      logical :: used
      integer :: i,j,k
      real,dimension(isd:ied,jsd:jed,nk_cbed)    :: cbed_tmask
      real,dimension(isd:ied,jsd:jed,nk_cbed+1)    :: cbed_tmask_i
      ! Make a cbed mask. Note: it seems grid_tmask(:,:,k) does not depend on k
      !do k=1,nk_cbed ; cbed_tmask(:,:,k) = grid_tmask(:,:,nk) ; enddo
      do j = jsc, jec; do i = isc, iec; do k=1,nk_cbed ;
               cbed_tmask(i,j,k) = grid_tmask(i,j,nk) ; enddo; enddo; enddo
      do j = jsc, jec; do i = isc, iec; do k=1,nk_cbed+1 ;
               cbed_tmask_i(i,j,k) = grid_tmask(i,j,nk) ; enddo; enddo; enddo


      used = send_data(cbed%id_o2, cbed%f_o2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_om1, cbed%f_om1, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_om2, cbed%f_om2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_om3, cbed%f_om3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_nh4, cbed%f_nh4, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_no3, cbed%f_no3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_dic, cbed%f_dic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_odu, cbed%f_odu, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_talk, cbed%f_talk, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_calc, cbed%f_calc, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_arag, cbed%f_arag, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_ca2, cbed%f_ca2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_po4, cbed%f_po4, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)

      ! diags
      ! 3D diags
      used = send_data(cbed%id_TOC, cbed%TOC, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_TIC, cbed%TIC, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_o2, cbed%R_om_o2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_no3, cbed%R_om_no3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_anaerobic, cbed%R_om_anaerobic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_dic, cbed%R_dic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_nox, cbed%R_nox, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_anammox, cbed%R_anammox, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_oduox, cbed%R_oduox, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_talk_org, cbed%R_talk_org, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_talk_inorg, cbed%R_talk_inorg, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_talk, cbed%R_talk, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_cbed_bioirri, cbed%cbed_bioirri, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_cbed_omega_calc, cbed%cbed_omega_calc, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_cbed_omega_arag, cbed%cbed_omega_arag, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_cbed_ph, cbed%cbed_ph, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      ! 2D diags
      used = send_data(cbed%id_o2_flux, cbed%o2_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_nh4_flux, cbed%nh4_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_no3_flux, cbed%no3_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_dic_flux, cbed%dic_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_talk_flux, cbed%talk_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_odu_flux, cbed%odu_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_burial_om, cbed%burial_om, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_denit, cbed%denit, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_k1, cbed%cbed_k1, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_k2, cbed%cbed_k2, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_k3, cbed%cbed_k3, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_w, cbed%cbed_w, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

      used = send_data(cbed%id_cbed_b_o2_acc, cbed%cbed_b_o2_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_dic_acc, cbed%cbed_b_dic_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_nh4_acc, cbed%cbed_b_nh4_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_no3_acc, cbed%cbed_b_no3_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_alk_acc, cbed%cbed_b_alk_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_odu_acc, cbed%cbed_b_odu_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_ca2_acc, cbed%cbed_b_ca2_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_b_po4_acc, cbed%cbed_b_po4_acc, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

      ! 1D diags
      used = send_data(cbed%id_dz_cbed, cbed%dz_cbed, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_z_cbed_mid, cbed%z_cbed_mid, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)

      ! 3D diags, CBED grid interfaces
      used = send_data(cbed%id_cbed_Db, cbed%cbed_Db, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_o2, cbed%cbed_D_o2, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in= jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_dic, cbed%cbed_D_dic, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_nh4, cbed%cbed_D_nh4, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_no3, cbed%cbed_D_no3, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_odu, cbed%cbed_D_odu, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_por, cbed%cbed_por, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_svf, cbed%cbed_svf, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)



   end subroutine generic_CBED_send_diagnostics

   subroutine generic_CBED_end()
      !Locals
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file
      type(domain2D),pointer :: domain

      !Resgister restarts
      call g_tracer_get_domain(domain)
      restart_file = 'RESTART/generic_CBED.res.nc'
      file_open_success=open_file(fileobj, trim(restart_file),"overwrite", domain, is_restart=.true.)
      if (file_open_success) then
         call register_axis(fileobj,'x','x')
         call register_axis(fileobj,'y','y')
         !< Register the domain decomposed dimensions as variables so that the combiner can work correctly
         call register_field(fileobj, "x", "double", (/"x"/))
         call register_field(fileobj, "y", "double", (/"y"/))
         call register_axis(fileobj,'lev',nk_cbed)
         ! register the restart variables
         call register_restart_field(fileobj, "cbed_o2", cbed%f_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om1", cbed%f_om1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om2", cbed%f_om2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om3", cbed%f_om3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_nh4", cbed%f_nh4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_no3", cbed%f_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_dic", cbed%f_dic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_odu", cbed%f_odu, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_talk", cbed%f_talk, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_calc", cbed%f_calc, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_arag", cbed%f_arag, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_ca2", cbed%f_ca2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_po4", cbed%f_po4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_b_o2_acc", cbed%cbed_b_o2_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_dic_acc", cbed%cbed_b_dic_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_nh4_acc", cbed%cbed_b_nh4_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_no3_acc", cbed%cbed_b_no3_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_alk_acc", cbed%cbed_b_alk_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_odu_acc", cbed%cbed_b_odu_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_ca2_acc", cbed%cbed_b_ca2_acc, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_b_po4_acc", cbed%cbed_b_po4_acc, (/"x","y"/))

         call write_restart(fileobj)
         call close_file(fileobj)
      else
         call error_mesg( 'generic_CBED_end', 'Cannot open restarts for write.', FATAL )
      endif

      !Deallocate arrays
      deallocate(cbed%f_o2)
      deallocate(cbed%f_om1)
      deallocate(cbed%f_om2)
      deallocate(cbed%f_om3)
      deallocate(cbed%f_nh4)
      deallocate(cbed%f_no3)
      deallocate(cbed%f_dic)
      deallocate(cbed%f_odu)
      deallocate(cbed%f_talk)
      deallocate(cbed%f_calc)
      deallocate(cbed%f_arag)
      deallocate(cbed%f_ca2)
      deallocate(cbed%f_po4)

      ! diags
      ! 3D diags
      deallocate(cbed%TOC)
      deallocate(cbed%TIC)
      deallocate(cbed%R_om_o2)
      deallocate(cbed%R_om_no3)
      deallocate(cbed%R_om_anaerobic)
      deallocate(cbed%R_dic)
      deallocate(cbed%R_nox)
      deallocate(cbed%R_anammox)
      deallocate(cbed%R_oduox)
      deallocate(cbed%R_talk_org)
      deallocate(cbed%R_talk_inorg)
      deallocate(cbed%R_talk)
      deallocate(cbed%cbed_bioirri)
      deallocate(cbed%cbed_omega_calc)
      deallocate(cbed%cbed_omega_arag)
      deallocate(cbed%cbed_ph)
      !2D diags
      deallocate(cbed%o2_flux)
      deallocate(cbed%nh4_flux)
      deallocate(cbed%no3_flux)
      deallocate(cbed%dic_flux)
      deallocate(cbed%talk_flux)
      deallocate(cbed%odu_flux)
      deallocate(cbed%burial_om)
      deallocate(cbed%denit)
      deallocate(cbed%cbed_k1)
      deallocate(cbed%cbed_k2)
      deallocate(cbed%cbed_k3)
      deallocate(cbed%cbed_w)
      deallocate(cbed%cbed_b_o2_acc)
      deallocate(cbed%cbed_b_dic_acc)
      deallocate(cbed%cbed_b_nh4_acc)
      deallocate(cbed%cbed_b_no3_acc)
      deallocate(cbed%cbed_b_alk_acc)
      deallocate(cbed%cbed_b_odu_acc)
      deallocate(cbed%cbed_b_ca2_acc)
      deallocate(cbed%cbed_b_po4_acc)


      ! 3D diags, CBED grid
      deallocate(cbed%dz_cbed)
      deallocate(cbed%z_cbed_mid)

      deallocate(por)
      deallocate(svf)

      deallocate(w)
      deallocate(Db_0)
      deallocate(Db)
      deallocate(bioirri_0)
      deallocate(bioirri)

      deallocate(D_o2)
      deallocate(D_dic)
      deallocate(D_nh4)
      deallocate(D_no3)
      deallocate(D_odu)
      deallocate(D_ca2)
      deallocate(D_po4)

      deallocate(k1)
      deallocate(k2)
      deallocate(k3)

      ! 3D diags, CBED grid interfaces
      deallocate(cbed%cbed_Db)
      deallocate(cbed%cbed_D_o2)
      deallocate(cbed%cbed_D_dic)
      deallocate(cbed%cbed_D_nh4)
      deallocate(cbed%cbed_D_no3)
      deallocate(cbed%cbed_D_odu)
      deallocate(cbed%cbed_por)
      deallocate(cbed%cbed_svf)


   end subroutine generic_CBED_end

   !-----------------------------------
   !! Tridiagonal matrix coefficients for tridiagonal solver.
   ! a(k)= (-ea(i,j,k)-sink(i,j,k))/h_old(k)
   ! b(k)= (h_old(k)+eb(i,j,k)+ea(i,j,k)+sink(i,j,k+1))/h_old(k)
   ! c(k)= -eb(i,j,k)/h_old(k)
   !----------------------------------------

   !! subroutine vertdiff_CBED !!
   ! VF issue is fixed. logic verified in R, and compared with R model.
   ! The solver logic reproduces the results of a diagenetic model solved through ReacTran package in R.

   subroutine vertdiff_CBED(cobalt_tracer_list, cobalt, cbed_field, field_name, D, w, VF, grid_kmt, dt, tau, isc, iec, jsc, jec, isd, ied, jsd, jed, nk, nk_cbed)
      ! Arguments
      type(g_tracer_type),          pointer       :: cobalt_tracer_list
      type(generic_COBALT_type),    intent(inout) :: cobalt
      real, dimension(:,:,:),       intent(inout) :: cbed_field  ! cbed tracer concentration field
      character(len=*),             intent(in)    :: field_name  ! Name of the cbed field
      real, dimension(:,:,:),       intent(in)    :: D           ! diffusion
      real, dimension(:,:,:),       intent(in)    :: w           ! sinking velocity or sedimentation rate
      real, dimension(:,:,:),       intent(in)    :: VF          ! volume fraction
      integer, dimension(:,:),      intent(in)    :: grid_kmt
      real,                         intent(in)    :: dt
      integer,                      intent(in)    :: tau
      integer,                      intent(in)    :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed

      ! Locals
      integer :: i, j, k
      logical :: is_solid

      ! Tridiagonal matrix arrays (1D column)
      real, dimension(nk_cbed)   :: a, b, c, f_old

      ! Transport intermediate arrays (1D column)
      real, dimension(nk_cbed+1) :: dist          ! True distance between nodes
      real, dimension(nk_cbed)   :: VF_cell       ! Cell-centered volume fraction
      real, dimension(nk_cbed)   :: capacity      ! Phase volume of the cell
      real, dimension(nk_cbed+1) :: K_diff        ! Diffusive conductance
      real, dimension(nk_cbed+1) :: K_adv         ! Advective conductance

      real :: DiffIn, DiffOut, AdvIn, AdvOut
      real :: btm_tracer_conc

      ! Local parameters. Fractions of total organic matter flux assigned to each reactivity class.
      real, dimension(isc:iec,jsc:jec) :: frac_OM1, frac_OM2, frac_OM3
      !real, parameter :: frac_OM1 = 0.70
      !real, parameter :: frac_OM2 = 0.20
      !real, parameter :: frac_OM3 = 0.10

      ! -----------------------------------------------------------------------
      ! 1. Determine if the tracer is a solid or a solute
      ! -----------------------------------------------------------------------
      is_solid = .false.
      if (trim(field_name) == "f_om1" .or. &
         trim(field_name) == "f_om2" .or. &
         trim(field_name) == "f_om3" .or. &
         trim(field_name) == "f_calc" .or. &
         trim(field_name) == "f_arag") then
         is_solid = .true.
      endif

      ! -----------------------------------------------------------------------
      ! 2. Compute true distances between cell centers for gradients
      ! -----------------------------------------------------------------------
      dist(1) = dz_cbed(1) / 2.0
      do k = 2, nk_cbed
         dist(k) = 0.5 * (dz_cbed(k-1) + dz_cbed(k))
      enddo
      dist(nk_cbed+1) = dz_cbed(nk_cbed) / 2.0

      ! -----------------------------------------------------------------------
      ! 3. Main spatial loops
      ! -----------------------------------------------------------------------
      ! calculate fractions of total organic matter flux assigned to each reactivity class as a func of bathymetric depth (m)
      if (is_solid) then
         do j = jsc, jec; do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  if (cbed%use_depth_dependent_OM_frac) then
                     frac_OM1(i,j) = min(0.80, 0.65*(100.0/cobalt%zt(i,j,nk))**0.5)
                     frac_OM3(i,j) = max(0.03, 0.04*(100.0/cobalt%zt(i,j,nk))**(-0.3))
                     frac_OM2(i,j) = 1.0 - (frac_OM1(i,j) + frac_OM3(i,j))
                  else
                     frac_OM1(i,j) = 0.70
                     frac_OM2(i,j) = 0.20
                     frac_OM3(i,j) = 0.10
                  endif
               endif
            enddo; enddo
      endif

      do j = jsc, jec
         do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then

               ! --- A. Get bottom water concentration ---
               btm_tracer_conc = 0.0
               if (trim(field_name) == "f_o2") then
                  btm_tracer_conc = max(0.0, cobalt%btm_o2(i,j) * cobalt%Rho_0)
               else if (trim(field_name) == "f_nh4") then
                  btm_tracer_conc = max(0.0, cobalt%f_nh4(i,j,nk) * cobalt%Rho_0)
               else if (trim(field_name) == "f_no3") then
                  btm_tracer_conc = max(0.0, cobalt%btm_no3(i,j) * cobalt%Rho_0)
               else if (trim(field_name) == "f_dic") then
                  btm_tracer_conc = max(0.0, cobalt%btm_dic(i,j) * cobalt%Rho_0)
               else if (trim(field_name) == "f_odu") then
                  btm_tracer_conc = abs(min(0.0, cobalt%btm_o2(i,j) * cobalt%Rho_0))
               else if (trim(field_name) == "f_talk") then
                  btm_tracer_conc = max(0.0, cobalt%btm_alk(i,j) * cobalt%Rho_0)
               else if (trim(field_name) == "f_ca2") then
                  btm_tracer_conc = max(0.0, fn_sw_ca2(cobalt%btm_salt(i,j), cobalt%Rho_0))
               else if (trim(field_name) == "f_po4") then
                  btm_tracer_conc = max(0.0, cobalt%f_po4(i,j,nk) * cobalt%Rho_0)
               endif

               ! --- B. Compute Cell Capacities and Interface Conductances ---
               do k = 1, nk_cbed
                  ! Average volume fraction to the cell center
                  VF_cell(k) = 0.5 * (VF(i,j,k) + VF(i,j,k+1))
                  ! Phase volume capacity
                  capacity(k) = VF_cell(k) * dz_cbed(k)
               enddo

               do k = 1, nk_cbed+1
                  ! Diffusive conductance (D * VF / distance)
                  K_diff(k) = VF(i,j,k) * D(i,j,k) / dist(k)
                  ! Advective conductance (Upwind: w * VF)
                  K_adv(k)  = max(0.0, VF(i,j,k) * w(i,j,k))
               enddo

               ! --- C. Build Tridiagonal Matrix ---
               do k = 1, nk_cbed

                  ! Transport rates normalized by cell capacity
                  DiffIn  = K_diff(k)   * dt / capacity(k)
                  DiffOut = K_diff(k+1) * dt / capacity(k)
                  AdvIn   = K_adv(k)    * dt / capacity(k)
                  AdvOut  = K_adv(k+1)  * dt / capacity(k)

                  if (k == 1) then
                     ! Top Boundary
                     a(1) = 0.0
                     c(1) = -DiffOut

                     if (.not. is_solid) then
                        ! Solutes: Robin BC
                        b(1) = 1.0 + DiffIn + DiffOut + AdvOut
                        f_old(1) = cbed_field(i,j,1) + (DiffIn + AdvIn) * btm_tracer_conc
                     else
                        ! Solids: Particle rain
                        b(1) = 1.0 + DiffOut + AdvOut
                        f_old(1) = cbed_field(i,j,1)

                        ! Add organic matter/solid fluxes normalized by capacity
                        if (trim(field_name) == "f_om1") then
                           f_old(1) = f_old(1) + (frac_OM1(i,j) * cobalt%fntot_btm(i,j) * cobalt%c_2_n * dt) / capacity(1)
                        else if (trim(field_name) == "f_om2") then
                           f_old(1) = f_old(1) + (frac_OM2(i,j) * cobalt%fntot_btm(i,j) * cobalt%c_2_n * dt) / capacity(1)
                        else if (trim(field_name) == "f_om3") then
                           f_old(1) = f_old(1) + (frac_OM3(i,j) * cobalt%fntot_btm(i,j) * cobalt%c_2_n * dt) / capacity(1)
                        else if (trim(field_name) == "f_calc") then
                           f_old(1) = f_old(1) + (cobalt%f_cadet_calc_btf(i,j,1) * dt) / capacity(1)
                        else if (trim(field_name) == "f_arag") then
                           f_old(1) = f_old(1) + (cobalt%f_cadet_arag_btf(i,j,1) * dt) / capacity(1)
                        endif
                     endif

                  else if (k == nk_cbed) then
                     ! Bottom Boundary
                     a(k) = -(DiffIn + AdvIn)
                     b(k) = 1.0 + DiffIn + AdvOut
                     c(k) = 0.0
                     f_old(k) = cbed_field(i,j,k)

                  else
                     ! Interior Nodes
                     a(k) = -(DiffIn + AdvIn)
                     b(k) = 1.0 + DiffIn + DiffOut + AdvOut
                     c(k) = -DiffOut
                     f_old(k) = cbed_field(i,j,k)
                  endif
               enddo

               ! --- D. Solve the System ---
               call CBED_tridag_solver_Press_et_al(a, b, c, f_old, cbed_field(i,j,:), nk_cbed)

            endif
         enddo
      enddo

   end subroutine vertdiff_CBED

!!!! Copied from Niki's code.
   subroutine CBED_tridag_solver_Press_et_al(a,b,c,r,u,n)
      integer, intent(in) :: n
      real,    intent(in) :: a(n),b(n),c(n),r(n)
      real,    intent(inout) :: u(n)
      real    :: bet,gam(n)
      integer :: k
      bet=b(1)
      u(1)=r(1)/bet
      do k=2,n
         gam(k)=c(k-1)/bet
         bet=b(k)-a(k)*gam(k)
         u(k)=(r(k)-a(k)*u(k-1))/bet
      enddo
      do k=n-1,1,-1
         u(k)=u(k)-gam(k+1)*u(k+1)
      enddo
   end subroutine CBED_tridag_solver_Press_et_al



   subroutine generic_CBED_update_from_source(cobalt_tracer_list, cobalt, phyto, ilb, jlb, mask_coast, &
      grid_tmask, grid_kmt, isc,iec, jsc,jec, isd,ied, jsd,jed, nk, r_dt, dt, tau, model_time, rho_dzt, dzt, internal_heat)

      type(g_tracer_type),          pointer       :: cobalt_tracer_list
      type(generic_COBALT_type),    intent(inout) :: cobalt
      type(phytoplankton), dimension(NUM_PHYTO), intent(inout) :: phyto
      integer,                      intent(in)    :: ilb, jlb
      real, dimension(:,:,:),       intent(in)    :: grid_tmask
      integer, dimension(:,:),      intent(in)    :: mask_coast, grid_kmt
      integer,                      intent(in)    :: isc,iec, jsc,jec, isd,ied, jsd,jed, nk
      real,                         intent(in)    :: r_dt, dt
      integer,                      intent(in)    :: tau
      type(time_type),              intent(in)    :: model_time
      real, dimension(ilb:,jlb:,:), intent(in)    :: rho_dzt, dzt
      real, dimension(ilb:,jlb:),   intent(in), optional :: internal_heat

      integer :: i, j, k
      integer :: stdoutunit
      real :: fpoc_btm, drho_dzt, log10_fpoc_btm
      integer, dimension(isc:iec,jsc:jec) :: k_bot
      real,    dimension(isc:iec,jsc:jec) :: rho_dzt_bot

      ! local parameters for bgc reactions
      ! Local parameters. Fractions of total organic matter flux assigned to each reactivity class. (only needed in vertdiff CBED, writing here for testing)
      real, dimension(isc:iec,jsc:jec) :: frac_OM1, frac_OM2, frac_OM3
      !real, parameter :: frac_OM1 = 0.70
      !real, parameter :: frac_OM2 = 0.20
      !real, parameter :: frac_OM3 = 0.10

      real, dimension(isc:iec,jsc:jec) :: cbed_burial_frac  ! = (OM burial / OM rain)  ! replaces cobalt%burial_frac
      real, dimension(isc:iec,jsc:jec) :: cbed_org_alk  ! organic alkalinity production from OM degradation


      real, parameter :: k_adj_denit = 0.1
      real, parameter :: k_adj_anoxia = 0.001

      real, parameter :: ks_o2 = 0.008   ! O2 half saturation constant (mol/m3)
      real, parameter :: ks_no3 = 0.001  ! NO3 half saturation constant (mol/m3)

      real, parameter :: k_nox = (2.0*10.0**5.0)/spery   ! 2e5 ! mol-1 m3 s-1 (from the original: mmol-1 L yr-1) !nitrification rate constant
      real, parameter :: k_ana = (10.0**5.0) /spery     !   ! 1e5               !anammox rate constant
      real, parameter :: k_oduox = (10.0**5.0) /spery   !              1e6      !ODU oxidation rate constant

      real, parameter :: Q10 = 1.88
      real, dimension(isc:iec,jsc:jec) :: Q10_factor

      real, dimension(isc:iec,jsc:jec) :: p_2_c  ! P:C ratio is calculated, because P:C or P:N is not in Redfield ratio in COBALT

      ! Reaction rates
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_o2, R_om2_o2, R_om3_o2
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_no3, R_om2_no3, R_om3_no3
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_anoxic, R_om2_anoxic, R_om3_anoxic
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_dic_om1, R_dic_om2, R_dic_om3
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_nox, R_ana, R_oduox, odu_depo
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_talk_org, R_talk_inorg, R_talk

      ! local variables for positive concentration constraint and for calculating reaction rates based on that.
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: c_om1, c_om2, c_om3, c_o2, c_no3, c_nh4, c_dic, c_odu, c_talk, c_calc, c_arag, c_ca2, c_po4

      ! b terms
      real, dimension(isc:iec,jsc:jec) :: b_o2, b_dic, b_nh4, b_no3, b_odu, b_alk, b_ca2, b_po4
      real, dimension(isc:iec,jsc:jec) :: b_o2_sub, b_dic_sub, b_nh4_sub, b_no3_sub, b_odu_sub, b_alk_sub, b_ca2_sub, b_po4_sub

      ! variables for sub-stepping reactions | adaptive time stepping for reactions
      integer :: n_req_o2, n_req_no3, n_req_nh4, n_req_odu
      !real    :: dt_sub
      real    :: max_o2_sink, max_no3_sink, max_nh4_sink, max_odu_sink
      integer :: sub_step
      integer, dimension(isc:iec,jsc:jec) :: n_sub
      real, dimension(isc:iec,jsc:jec) :: dt_sub

      ! carbonate system related variables
      real :: omega_arag_crit, omega_calc_crit
      real :: n_diss_arag_gt_crit, n_diss_arag_lt_crit, n_diss_calc_gt_crit, n_diss_calc_lt_crit
      real :: n_diss_calc, n_diss_arag, n_prec_calc
      real :: k_diss_arag_gt_crit, k_diss_arag_lt_crit, k_diss_calc_gt_crit, k_diss_calc_lt_crit 
      real :: k_diss_calc, k_diss_arag, k_prec_calc, k_prec_arag

      real, dimension(isc:iec,jsc:jec,nk_cbed) :: cbed_omega_arag, cbed_omega_calc
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_diss_calc, R_diss_arag, R_prec_calc, R_prec_arag
      real :: c_dic_co2, c_po4_co2, c_talk_co2, c_nh4_co2, c_ca2_co2, c_h2s_co2, c_sio4_co2, htotal_dummy



      ! write the grid layers to register as diags
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed
                  cbed%dz_cbed(i,j,k) = dz_cbed(k)
                  cbed%z_cbed_mid(i,j,k) = z_cbed_mid(k)
               enddo
            endif
         enddo; enddo

      ! Porosity and solid volume fraction read from file.
      if (cbed%read_porosity_from_file) then
         call data_override('OCN', 'por', por(isc:iec,jsc:jec,1), model_time)

         do j = jsc, jec
            do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  do k = 2, nk_cbed+1
                     por(i,j,k) = por(i,j,1)
                  enddo

                  do k = 1, nk_cbed+1
                     svf(i,j,k) = 1.0 - por(i,j,k)
                  enddo

               endif
            enddo
         enddo
      else
         do j = jsc, jec
            do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  do k = 1, nk_cbed+1
                     por(i,j,k) = 0.8
                     svf(i,j,k) = 0.2
                  enddo
               endif
            enddo
         enddo
         !por = 0.8
         !svf = 0.2
      endif


      ! calculate required parameters in i,j loop
      do j = jsc, jec
         do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then

               !-------------------------
               ! calculate Q10 factor
               !-------------------------
               Q10_factor(i,j) = Q10**( (cobalt%btm_temp(i,j)-4.0)/10.0 )

               !-------------------------
               ! calculate P:C ratio
               !-------------------------
               p_2_c(i,j) = cobalt%fptot_btm(i,j)/(cobalt%fntot_btm(i,j)*cobalt%c_2_n + epsln)

               !----------------------------------
               ! Sedimentation rate calculation
               !----------------------------------
               do k = 1, nk_cbed+1

                  ! w (sedimentation rate, cm/year) m/s
                  w(i,j,k) = ( (cobalt%f_cadet_arag_btf(i,j,1)*100.0/2.71 + &
                     cobalt%f_cadet_calc_btf(i,j,1)*100.0/2.94 + &
                     cobalt%fsitot_btm(i,j)*60.0/2.65 + &
                     cobalt%f_lithdet_btf(i,j,1)/2.65 + &
                     cobalt%ffetot_btm(i,j)*160.0/5.24 + &
                     cobalt%fptot_btm(i,j)*120.0/2.3 + &
                     cobalt%fntot_btm(i,j)*cobalt%c_2_n*22.4/0.9)/10000.0*spery/svf(i,j,k) )/100.0/spery
               enddo

               !------------------
               !Bioturbation
               !------------------
               ! relation from Archer. POC flux unit in umol cm-2 y-1.
               Db_0(i,j) = ( 0.0232*((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**0.85) ) /1e4/spery ! in cobalt unit m2/s

               do k = 1, nk_cbed+1
                  ! relation from Archer. POC flux unit in umol cm-2 y-1.
                  Db(i,j,k) = max(0.0, Db_0(i,j)*exp(-(z_cbed_int(k)/Db_l)**2.0)*(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)/(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)+(20.0/1e3))) )
               enddo

               !--------------------
               ! Bioirrigation
               !--------------------
               ! relation from Archer. POC flux unit in umol cm-2 y-1.
               bioirri_0(i,j) = ( 11.0*(((atan((5.0*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery) -400.0)/400.0))/pi)+0.5) &
                  - 0.9 + 20.0*((max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0))/(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)+0.01)) * exp(-max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)/0.01) * &
                  ((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)/((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)+30.0)) )/spery   ! in cobalt unit s^-1
               do k = 1, nk_cbed
                  ! relation from Archer. POC flux unit in umol cm-2 y-1.
                  bioirri(i,j,k) = max(0.0, bioirri_0(i,j)*exp(-(z_cbed_mid(k)/bioirri_l)**2.0) )
               enddo

               !-----------------------------------
               !Calculate diffusion coefficients
               !-----------------------------------
               do k = 1, nk_cbed+1
                  D_o2(i,j,k)  = ( (0.031558+0.001428*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)   ! m2/s
                  D_dic(i,j,k) = ( (0.015179+0.000795*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_nh4(i,j,k) = ( (0.030926+0.001225*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_no3(i,j,k) = ( (0.030863+0.001153*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_odu(i,j,k) = ( (0.028938+0.001314*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_ca2(i,j,k) = ( (0.011771+0.000529*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_po4(i,j,k) = ( (0.009783+0.000513*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)

                  ! b_odu(i,j) = por(i,j,1)*D_odu(i,j,1)*((0.0 - max(0.0,cbed%f_odu(i,j,1)))/(dz_cbed(1)/2.0)) + &
                  ! sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(0.0 - max(0.0,c_odu(i,j,:))) )

                  ! if ( b_odu(i,j) .lt. -1.0e6 ) then
                  !    D_odu(i,j,k) = 0.01 * D_odu(i,j,k)
                  ! end if


               enddo

               !------------------
               !Calculate OM decay rates. k1,k2,k3 [unit: s-1]
               !------------------
               ! POC flux unit in umol cm-2 y-1. Unit of k is y-1
               k1(i,j) = ( 0.15*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**(0.85) )/spery
               k2(i,j) = ( 0.0015*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**(0.85) )/spery
               k3(i,j) = ( 0.00009*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**(0.85) )/spery


               !-----------------------
               ! Set concentrations to be positive to avoid negative reaction rates.
               ! This can happen when the tracer is very low and the change is large in one time step,
               ! causing it to go negative in the next time step. Setting to zero prevents runaway reactions from negative concentrations.
               ! Setting temporary variables to do that without changing the actual concentration fields, which are used in the diffusion-advection solver.
               ! This way we can prevent negative concentrations from causing unrealistic reaction rates, while still allowing the diffusion-advection solver to bring the concentration back up in the next time step.
               !------------------------
               do k = 1, nk_cbed
                  c_om1(i,j,k) = max(0.0, cbed%f_om1(i,j,k))
                  c_om2(i,j,k) = max(0.0, cbed%f_om2(i,j,k))
                  c_om3(i,j,k) = max(0.0, cbed%f_om3(i,j,k))
                  c_o2(i,j,k)  = max(0.0, cbed%f_o2(i,j,k))
                  c_no3(i,j,k) = max(0.0, cbed%f_no3(i,j,k))
                  c_nh4(i,j,k) = max(0.0, cbed%f_nh4(i,j,k))
                  c_dic(i,j,k) = max(0.0, cbed%f_dic(i,j,k))
                  c_odu(i,j,k) = max(0.0, cbed%f_odu(i,j,k))
                  c_talk(i,j,k) = max(0.0, cbed%f_talk(i,j,k))
                  c_calc(i,j,k) = max(0.0, cbed%f_calc(i,j,k))
                  c_arag(i,j,k) = max(0.0, cbed%f_arag(i,j,k))
                  c_ca2(i,j,k) = max(0.0, cbed%f_ca2(i,j,k))
                  c_po4(i,j,k) = max(0.0, cbed%f_po4(i,j,k))
               enddo

               !---------------------
               ! calculate reaction rates
               !----------------------
               do k = 1, nk_cbed

                  ! O₂ reaction rates
                  R_om1_o2(i,j,k) = k1(i,j)*c_om1(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  R_om2_o2(i,j,k) = k2(i,j)*c_om2(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  R_om3_o2(i,j,k) = k3(i,j)*c_om3(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  ! NO₃ reaction rates
                  R_om1_no3(i,j,k) = k_adj_denit*k1(i,j)*c_om1(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  R_om2_no3(i,j,k) = k_adj_denit*k2(i,j)*c_om2(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  R_om3_no3(i,j,k) = k_adj_denit*k3(i,j)*c_om3(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  ! ODU reaction rates
                  R_om1_anoxic(i,j,k) = k_adj_anoxia*k1(i,j)*c_om1(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  R_om2_anoxic(i,j,k) = k_adj_anoxia*k2(i,j)*c_om2(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                  R_om3_anoxic(i,j,k) = k_adj_anoxia*k3(i,j)*c_om3(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)

                  ! dic
                  R_dic_om1(i,j,k) = (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k))
                  R_dic_om2(i,j,k) = (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k))
                  R_dic_om3(i,j,k) = (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k))
                  ! nitrification
                  R_nox(i,j,k) = k_nox*c_nh4(i,j,k)*c_o2(i,j,k) * Q10_factor(i,j)
                  ! anammox
                  R_ana(i,j,k) = k_ana*c_nh4(i,j,k)*c_no3(i,j,k) * Q10_factor(i,j) !* (ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                  ! ODU oxidation
                  R_oduox(i,j,k) = k_oduox*c_odu(i,j,k)*c_o2(i,j,k) * Q10_factor(i,j)
                  odu_depo(i,j,k) = (R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))*min(1.0, 0.233*(w(i,j,k)*100.0*spery)**0.336)

                  ! TA calculation
                  R_talk_org(i,j,k) = svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n - p_2_c(i,j))*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) + &
                     svf(i,j,k)/por(i,j,k)*(0.8+1.0/cobalt%c_2_n - p_2_c(i,j))*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
                     svf(i,j,k)/por(i,j,k)*(1.0+1.0/cobalt%c_2_n - p_2_c(i,j))*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
                     2.0*R_nox(i,j,k) - 1.0*R_oduox(i,j,k) - 0.4*R_ana(i,j,k)

                  R_talk_inorg(i,j,k) = 0.0

                  R_talk(i,j,k) = R_talk_org(i,j,k) + R_talk_inorg(i,j,k)

                  ! calculations for diagnostics
                  cbed%R_om_o2(i,j,k) = R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)
                  cbed%R_om_no3(i,j,k) = R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)
                  cbed%R_om_anaerobic(i,j,k) = R_om1_anoxic(i,j,k) + R_om2_anoxic(i,j,k) + R_om3_anoxic(i,j,k)
                  cbed%R_dic(i,j,k) = R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)
                  cbed%R_nox(i,j,k) = R_nox(i,j,k)
                  cbed%R_anammox(i,j,k) = R_ana(i,j,k)
                  cbed%R_oduox(i,j,k) = R_oduox(i,j,k)
                  cbed%R_talk_org(i,j,k) = R_talk_org(i,j,k)
                  cbed%R_talk_inorg(i,j,k) = R_talk_inorg(i,j,k)
                  cbed%R_talk(i,j,k) = R_talk(i,j,k)

                  cbed%cbed_bioirri(i,j,k) = bioirri(i,j,k) ! bioirrigation diagnostics
               enddo


               !-------------------
               ! pass the b_terms calculated at t-1 timestep as b_terms
               !---------------------
               if (cbed%do_adaptive_time_stepping) then

                  b_o2(i,j) = cbed%cbed_b_o2_acc(i,j)
                  b_dic(i,j) = cbed%cbed_b_dic_acc(i,j)
                  b_nh4(i,j) = cbed%cbed_b_nh4_acc(i,j)
                  b_no3(i,j) = cbed%cbed_b_no3_acc(i,j)
                  b_odu(i,j) = cbed%cbed_b_odu_acc(i,j)
                  b_alk(i,j) = cbed%cbed_b_alk_acc(i,j)

               else
                  b_o2(i,j) = por(i,j,1)*D_o2(i,j,1)*((max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,1))/(dz_cbed(1)/2.0)) + &
                     por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) + &
                     sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,:)) )

                  b_nh4(i,j) = por(i,j,1)*D_nh4(i,j,1)*((max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,1))/(dz_cbed(1)/2.0)) + &
                     por(i,j,1)*w(i,j,1)*max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) + &
                     sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,:)) )

                  b_no3(i,j) = por(i,j,1)*D_no3(i,j,1)*((max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,1))/(dz_cbed(1)/2.0)) + &
                     por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) + &
                     sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,:)) )

                  b_dic(i,j) = por(i,j,1)*D_dic(i,j,1)*((max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,1))/(dz_cbed(1)/2.0)) + &
                     por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) + &
                     sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,:)) )

                  b_odu(i,j) = por(i,j,1)*D_odu(i,j,1)*((abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,1))/(dz_cbed(1)/2.0)) + &
                     por(i,j,1)*w(i,j,1)*abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) + &
                     sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,:)) )

                  b_alk(i,j) = sum(dz_cbed(:)*por(i,j,1:nk_cbed)*R_talk(i,j,:))
               endif


               !---------
               ! save the benthic fluxes as diagnostics. 2D diag
               !----------
               cbed%o2_flux(i,j)  = b_o2(i,j)
               cbed%nh4_flux(i,j) = b_nh4(i,j)
               cbed%no3_flux(i,j) = b_no3(i,j)
               cbed%dic_flux(i,j) = b_dic(i,j)

               cbed%talk_flux(i,j) = b_alk(i,j)
               ! cbed%talk_flux(i,j) = por(i,j,1)*D_dic(i,j,1)*((max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) - c_talk(i,j,1))/(dz_cbed(1)/2.0)) + &
               !    por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) + &
               !    sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) - c_talk(i,j,:)) )

               cbed%odu_flux(i,j) = b_odu(i,j)


               !-----------------
               !TOC and TIC (total organic carbon, total inorganic (calc+arag) carbon) diagnostics
               !----------------
               do k = 1, nk_cbed
                  cbed%TOC(i,j,k)  = (c_om1(i,j,k)+c_om2(i,j,k)+c_om3(i,j,k))*12.0/1e6/rho_s*100.0
                  cbed%TIC(i,j,k)  = (c_calc(i,j,k)+c_arag(i,j,k))*100.0/1e6/rho_s*100.0
               enddo

               !-----------------
               ! some other diags , other local variables
               !-----------------
               cbed%burial_om(i,j)  = (c_om1(i,j,nk_cbed)+c_om2(i,j,nk_cbed)+c_om3(i,j,nk_cbed))*w(i,j,nk_cbed+1) * svf(i,j,nk_cbed+1) ! mol/m2/s

               cbed%denit(i,j) = sum(dz_cbed(:)*(svf(i,j,1:nk_cbed)*0.8*cbed%R_om_no3(i,j,:) + por(i,j,1:nk_cbed)*1.6*R_ana(i,j,:)))

               !cbed%cbed_k1(i,j) = k1(i,j)
               !cbed%cbed_k2(i,j) = k2(i,j)
               !cbed%cbed_k3(i,j) = k3(i,j)
               cbed%cbed_w(i,j) = w(i,j,1)

               !------------------------------
               ! some other local variables
               !------------------------------
               cbed_burial_frac(i,j) = max(0.0, cbed%burial_om(i,j)/(cobalt%fntot_btm(i,j)*cobalt%c_2_n + epsln) )  ! burial / rain . ratio

               cbed_org_alk(i,j) = sum(dz_cbed(:)*por(i,j,1:nk_cbed)*R_talk(i,j,:))  ! mol m-2 s-1 (net production of alklinity from organic matter degradation)


               !--------------------------
               ! some more diags
               !--------------------------
               do k = 1, nk_cbed+1
                  cbed%cbed_Db(i,j,k) = Db(i,j,k)
                  cbed%cbed_D_o2(i,j,k) = D_o2(i,j,k)
                  cbed%cbed_D_nh4(i,j,k) = D_nh4(i,j,k)
                  cbed%cbed_D_no3(i,j,k) = D_no3(i,j,k)
                  cbed%cbed_D_dic(i,j,k) = D_dic(i,j,k)
                  cbed%cbed_D_odu(i,j,k) = D_odu(i,j,k)
                  cbed%cbed_por(i,j,k) = por(i,j,k)
                  cbed%cbed_svf(i,j,k) = svf(i,j,k)
               enddo


            endif
         enddo
      enddo    ! i,j loop


      ! Calculate the "b terms" to feed into cobalt.
      ! The b terms are calulated based on the t-1 time step. This is to maintain mass balance
      ! with cobalt-cbed. The surface flux calculated in the vertdiff_CBED subroutine
      ! is the flux at t-1, so we need to use the t-1 concentration field to calculate the b term,
      ! which will update the bottom flux as the verdiff_G is not yet called for cobalt.

      ! The b terms (benthic fluxes at sediment-water interface) are calculated as the sum of diffusive flux, advective flux, and bioirrigation flux.
      ! The concentration gradient for the diffusive flux is calculated using the t-1 concentration field, and the bioirrigation flux is calculated using the t-1 concentration field and the bioirrigation rate at t-1.
      ! This way, we ensure that the b term represents the correct fluxes based on the t-1 state of the system.

      ! Negative values of b term would indicate fluxes from the sediment to the water column, while positive values would indicate fluxes from the water column to the sediment i.e. consumption by the sediment.

      ! The concentration fields are updated (to the t time step) once the vertdiff_CBED subroutine is called. (and for COBALT, when the vertdiff_G subroutine is called)

      ! further modification related to b_o2, b_dic etc is done when it is passed to cobalt%b_* at the end of the code.

      ! ! calculate the b_terms
      ! if ( .not. cbed%do_adaptive_time_stepping) then
      !    do j = jsc, jec; do i = isc, iec
      !          if (grid_kmt(i,j) .gt. 0) then

      !             b_o2(i,j) = por(i,j,1)*D_o2(i,j,1)*((max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,1))/(dz_cbed(1)/2.0)) + &
      !                por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) + &
      !                sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,:)) )

      !             b_nh4(i,j) = por(i,j,1)*D_nh4(i,j,1)*((max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,1))/(dz_cbed(1)/2.0)) + &
      !                por(i,j,1)*w(i,j,1)*max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) + &
      !                sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,:)) )

      !             b_no3(i,j) = por(i,j,1)*D_no3(i,j,1)*((max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,1))/(dz_cbed(1)/2.0)) + &
      !                por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) + &
      !                sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,:)) )

      !             b_dic(i,j) = por(i,j,1)*D_dic(i,j,1)*((max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,1))/(dz_cbed(1)/2.0)) + &
      !                por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) + &
      !                sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,:)) )

      !             b_odu(i,j) = por(i,j,1)*D_odu(i,j,1)*((abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,1))/(dz_cbed(1)/2.0)) + &
      !                por(i,j,1)*w(i,j,1)*abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) + &
      !                sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,:)) )

      !             b_alk_org(i,j) = sum(dz_cbed(:)*por(i,j,1:nk_cbed)*R_talk(i,j,:))

      !          endif
      !       enddo;enddo
      ! endif


      ! ! test for adaptive time-step off. delete this later
      ! !n_sub = 1
      ! cbed%cbed_b_o2_acc = 0.0
      ! cbed%cbed_b_dic_acc = 0.0
      ! cbed%cbed_b_nh4_acc = 0.0
      ! cbed%cbed_b_no3_acc = 0.0
      ! cbed%cbed_b_alk_org_acc = 0.0
      ! cbed%cbed_b_odu_acc = 0.0
      ! if ( .not. cbed%do_adaptive_time_stepping) then
      !    do j = jsc, jec
      !       do i = isc, iec
      !          if (grid_kmt(i,j) > 0) then
      !             cbed%cbed_b_o2_acc(i,j) = 0.0
      !             cbed%cbed_b_dic_acc(i,j) = 0.0
      !             cbed%cbed_b_nh4_acc(i,j) = 0.0
      !             cbed%cbed_b_no3_acc(i,j) = 0.0
      !             cbed%cbed_b_alk_org_acc(i,j) = 0.0
      !             cbed%cbed_b_odu_acc(i,j) = 0.0
      !          endif
      !       enddo
      !    enddo
      ! endif


      ! ask NIKI: if this would lead to different ans due to indexing and domain decomposition.
      if (cbed%do_adaptive_time_stepping) then

         ! --- ADAPTIVE TIME-STEPPING CALCULATION (O2, NO3, NH4, ODU) ---

         n_sub = 1 ! Default to 1 macro step

         do j = jsc, jec
            do i = isc, iec
               if (grid_kmt(i,j) > 0) then
                  do k = 1, nk_cbed

                     ! ---------------------------------------------------------
                     ! 1. OXYGEN CONSTRAINT
                     ! Sinks: Aerobic Respiration, Nitrification, ODU Oxidation
                     ! ---------------------------------------------------------
                     if (c_o2(i,j,k) > 1.0e-6) then
                        max_o2_sink = svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) + &
                           (2.0*R_nox(i,j,k)+R_oduox(i,j,k))

                        ! The maximum sink is the total amount of O2 that could be consumed in this time step based on the current concentrations and reaction rates.
                        ! We divide by 0.5*c_o2, to get the number of sub-steps needed to ensure that we don't drop the cell concentration more than half in one time step.
                        ! This is a conservative estimate to ensure we don't overshoot and get negative concentrations.
                        ! We then take the ceiling of this number to get the number of sub-steps needed to ensure that we don't consume more O2 than is available in any sub-step.
                        ! We do this for O2, NO3, and NH4 and take the maximum number of sub-steps required among the three to ensure that we don't violate any of the constraints.

                        n_req_o2 = ceiling( (max_o2_sink * dt) / (0.8 * c_o2(i,j,k)) )
                        if (n_req_o2 > n_sub(i,j)) n_sub(i,j) = n_req_o2
                     endif

                     ! ---------------------------------------------------------
                     ! 2. NITRATE CONSTRAINT
                     ! Sinks: Denitrification, Anammox
                     ! ---------------------------------------------------------
                     if (c_no3(i,j,k) > 1.0e-6) then
                        max_no3_sink = svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + R_ana(i,j,k)

                        n_req_no3 = ceiling( (max_no3_sink * dt) / (0.8 * c_no3(i,j,k)) )
                        if (n_req_no3 > n_sub(i,j)) n_sub(i,j) = n_req_no3
                     endif

                     ! ---------------------------------------------------------
                     ! 3. AMMONIUM CONSTRAINT
                     ! Sinks: Nitrification, Anammox
                     ! ---------------------------------------------------------
                     if (c_nh4(i,j,k) > 1.0e-6) then
                        max_nh4_sink = R_nox(i,j,k) + R_ana(i,j,k)

                        n_req_nh4 = ceiling( (max_nh4_sink * dt) / (0.8 * c_nh4(i,j,k)) )
                        if (n_req_nh4 > n_sub(i,j)) n_sub(i,j) = n_req_nh4
                     endif

                     ! ---------------------------------------------------------
                     ! 4. ODU CONSTRAINT
                     ! Sinks: ODU oxidation
                     ! ---------------------------------------------------------
                     if (c_odu(i,j,k) > 1.0e-6) then
                        max_odu_sink = R_oduox(i,j,k)

                        n_req_odu = ceiling( (max_odu_sink * dt) / (0.8 * c_odu(i,j,k)) )
                        if (n_req_odu > n_sub(i,j)) n_sub(i,j) = n_req_odu
                     endif

                  enddo
               endif
            enddo
         enddo

         ! Cap the maximum number of sub-steps to prevent the ESM from hanging
         ! Note: This cap can be increased depending on how aggressive the coastal fluxes may get.
         do j = jsc, jec
            do i = isc, iec
               n_sub(i,j) = min(n_sub(i,j), 10)
               dt_sub(i,j) = dt / real(n_sub(i,j))
            enddo
         enddo

         ! Initialize macro-step accumulators for benthic fluxes to the ocean
         cbed%cbed_b_o2_acc = 0.0
         cbed%cbed_b_dic_acc = 0.0
         cbed%cbed_b_nh4_acc = 0.0
         cbed%cbed_b_no3_acc = 0.0
         cbed%cbed_b_alk_acc = 0.0
         cbed%cbed_b_odu_acc = 0.0
         cbed%cbed_b_ca2_acc = 0.0
         cbed%cbed_b_po4_acc = 0.0


         ! --- BEGIN ADAPTIVE SUB-STEPPING LOOP ---
         !do sub_step = 1, n_sub

         ! 1. Re-evaluate positive concentrations for this specific sub-step
         do j = jsc, jec
            do i = isc, iec
               if (grid_kmt(i,j) > 0) then
                  do sub_step = 1, n_sub(i,j)
                     do k = 1, nk_cbed
                        c_om1(i,j,k) = max(0.0, cbed%f_om1(i,j,k))
                        c_om2(i,j,k) = max(0.0, cbed%f_om2(i,j,k))
                        c_om3(i,j,k) = max(0.0, cbed%f_om3(i,j,k))
                        c_o2(i,j,k)  = max(0.0, cbed%f_o2(i,j,k))
                        c_no3(i,j,k) = max(0.0, cbed%f_no3(i,j,k))
                        c_nh4(i,j,k) = max(0.0, cbed%f_nh4(i,j,k))
                        c_dic(i,j,k) = max(0.0, cbed%f_dic(i,j,k))
                        c_odu(i,j,k) = max(0.0, cbed%f_odu(i,j,k))
                        c_talk(i,j,k) = max(0.0, cbed%f_talk(i,j,k))
                        c_calc(i,j,k) = max(0.0, cbed%f_calc(i,j,k))
                        c_arag(i,j,k) = max(0.0, cbed%f_arag(i,j,k))
                        c_ca2(i,j,k) = max(0.0, cbed%f_ca2(i,j,k))
                        c_po4(i,j,k) = max(0.0, cbed%f_po4(i,j,k))


                        ! O₂ reaction rates
                        R_om1_o2(i,j,k) = k1(i,j)*c_om1(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        R_om2_o2(i,j,k) = k2(i,j)*c_om2(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        R_om3_o2(i,j,k) = k3(i,j)*c_om3(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        ! NO₃ reaction rates
                        R_om1_no3(i,j,k) = k_adj_denit*k1(i,j)*c_om1(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        R_om2_no3(i,j,k) = k_adj_denit*k2(i,j)*c_om2(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        R_om3_no3(i,j,k) = k_adj_denit*k3(i,j)*c_om3(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        ! ODU reaction rates
                        R_om1_anoxic(i,j,k) = k_adj_anoxia*k1(i,j)*c_om1(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        R_om2_anoxic(i,j,k) = k_adj_anoxia*k2(i,j)*c_om2(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
                        R_om3_anoxic(i,j,k) = k_adj_anoxia*k3(i,j)*c_om3(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)

                        ! dic
                        R_dic_om1(i,j,k) = (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k))
                        R_dic_om2(i,j,k) = (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k))
                        R_dic_om3(i,j,k) = (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k))
                        ! nitrification
                        R_nox(i,j,k) = k_nox*c_nh4(i,j,k)*c_o2(i,j,k) * Q10_factor(i,j)
                        ! anammox
                        R_ana(i,j,k) = k_ana*c_nh4(i,j,k)*c_no3(i,j,k) * Q10_factor(i,j) !* (ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                        ! ODU oxidation
                        R_oduox(i,j,k) = k_oduox*c_odu(i,j,k)*c_o2(i,j,k) * Q10_factor(i,j)
                        odu_depo(i,j,k) = (R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))*min(1.0, 0.233*(w(i,j,k)*100.0*spery)**0.336)


                        !------
                        ! start carbonate system calculations
                        !------
                        ! 1) call co2calc. get omega. 2) determine n_diss, k_diss etc. 3) write R_diss, R_prec


                        htotal_dummy = cobalt%f_htotal(i,j,nk)
                        
                        ! Convert mol/m3 -> mol/kg before supplying to FMS_co2calc. 
                        ! *_co2 are local variables used for this purpose. 
                        c_dic_co2 = c_dic(i,j,k) / cobalt%Rho_0
                        c_po4_co2 = c_po4(i,j,k) / cobalt%Rho_0
                        c_sio4_co2 = cobalt%f_sio4(i,j,nk)     ! comes from cobalt directly. already in mol/kg
                        c_talk_co2 = c_talk(i,j,k) / cobalt%Rho_0
                        c_nh4_co2 = c_nh4(i,j,k) / cobalt%Rho_0
                        c_h2s_co2 = 0.5*c_odu(i,j,k) / cobalt%Rho_0  !approximates H2S. 
                        c_ca2_co2 = c_ca2(i,j,k) / cobalt%Rho_0


                        call FMS_co2calc_point(mask=grid_tmask(i,j,nk),&
                           t_in=cobalt%btm_temp(i,j),                  &
                           s_in=cobalt%btm_salt(i,j),                  &
                           dic_in=c_dic_co2,                        &
                           pt_in=c_po4_co2,                          &
                           sit_in=c_sio4_co2,                         &
                           ta_in=c_talk_co2,                          &
                        !InOut
                           htotal=htotal_dummy,                       &
                        !Optional In
                           zt=cobalt%zt(i,j,nk),                          &
                           nh4_in=c_nh4_co2,                           &
                           h2s_in=c_h2s_co2,                       &
                           ca_in=c_ca2_co2,                            &
                           optCON_in='mol/kg',                            &
                        !OUT
                           omega_arag=cbed%cbed_omega_arag(i,j,k), &
                           omega_calc=cbed%cbed_omega_calc(i,j,k))

                        ! track cbed_htotal, convert to ph back outside subroutine. 
                        cbed%cbed_ph(i,j,k) = - LOG10(htotal_dummy)
                        


                        ! assign rate constants (RADI) {
                        omega_arag_crit = 0.835
                        omega_calc_crit = 0.828

                        n_diss_arag_gt_crit = 0.13  ! gt = greater than critical
                        n_diss_arag_lt_crit = 1.46  ! lt = less than critical

                        n_diss_calc_gt_crit = 0.11
                        n_diss_calc_lt_crit = 4.7

                        n_prec_calc = 1.76

                        k_diss_arag_gt_crit = (3.8*10.0**(-3.0))  !unit year-1
                        k_diss_arag_lt_crit = (4.2*10.0**(-2.0))

                        k_diss_calc_gt_crit = (6.3*10.0**(-3.0))
                        k_diss_calc_lt_crit = 20.0 

                        k_prec_calc = 0.4/spery  !unit mol m-3 year-1 converted to mol m-3 s-1
                        k_prec_arag = 0.0
                        
                        !--RADI--!}

                        ! assign the parameters according to critical saturation state
                        ! n_diss
                        if (cbed%cbed_omega_arag(i,j,k) >= omega_arag_crit .and. cbed%cbed_omega_arag(i,j,k) < 1.0) then
                           n_diss_arag = n_diss_arag_gt_crit
                        else
                           n_diss_arag = n_diss_arag_lt_crit
                        endif

                        if (cbed%cbed_omega_calc(i,j,k) >= omega_calc_crit .and. cbed%cbed_omega_calc(i,j,k) < 1.0) then
                           n_diss_calc = n_diss_calc_gt_crit
                        else
                           n_diss_calc = n_diss_calc_lt_crit
                        endif

                        ! k_diss
                        if (cbed%cbed_omega_arag(i,j,k) >= omega_arag_crit .and. cbed%cbed_omega_arag(i,j,k) < 1.0) then
                           k_diss_arag = k_diss_arag_gt_crit/spery
                        else
                           k_diss_arag = k_diss_arag_lt_crit/spery
                        endif

                        if (cbed%cbed_omega_calc(i,j,k) >= omega_calc_crit .and. cbed%cbed_omega_calc(i,j,k) < 1.0) then
                           k_diss_calc = k_diss_calc_gt_crit/spery
                        else
                           k_diss_calc = k_diss_calc_lt_crit/spery
                        endif


                        ! R_diss_calc, R_diss_arag, R_prec_calc, R_prec_arag

                        if (cbed%cbed_omega_arag(i,j,k) < 1.0) then
                           R_diss_arag(i,j,k) = Q10_factor(i,j) * k_diss_arag * c_arag(i,j,k) * (1.0-cbed%cbed_omega_arag(i,j,k))**n_diss_arag
                           R_prec_arag(i,j,k) = 0.0
                        else
                           R_diss_arag(i,j,k) = 0.0
                           R_prec_arag(i,j,k) = k_prec_arag * (cbed%cbed_omega_arag(i,j,k) - 1.0)
                        endif

                        if (cbed%cbed_omega_calc(i,j,k) < 1.0) then
                           R_diss_calc(i,j,k) = Q10_factor(i,j) * k_diss_calc * c_calc(i,j,k) * (1.0-cbed%cbed_omega_calc(i,j,k))**n_diss_calc
                           R_prec_calc(i,j,k) = 0.0
                        else
                           R_diss_calc(i,j,k) = 0.0
                           R_prec_calc(i,j,k) = k_prec_calc * (cbed%cbed_omega_calc(i,j,k) - 1.0)**n_prec_calc
                        endif

                        !
                        !end carbonate system calculations
                        !-----------


                        ! TA calculation
                        R_talk_org(i,j,k) = svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n - p_2_c(i,j))*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) + &
                           svf(i,j,k)/por(i,j,k)*(0.8+1.0/cobalt%c_2_n - p_2_c(i,j))*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
                           svf(i,j,k)/por(i,j,k)*(1.0+1.0/cobalt%c_2_n - p_2_c(i,j))*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
                           2.0*R_nox(i,j,k) - 1.0*R_oduox(i,j,k) - 0.4*R_ana(i,j,k)

                        R_talk_inorg(i,j,k) = 2.0*svf(i,j,k)/por(i,j,k)*(R_diss_arag(i,j,k) + R_diss_calc(i,j,k)) - &
                           2.0*(R_prec_arag(i,j,k) + R_prec_calc(i,j,k))

                        R_talk(i,j,k) = R_talk_org(i,j,k) + R_talk_inorg(i,j,k)

                     enddo !k

                     ! b terms for averaging throughout adaptive time stepping

                     b_o2_sub(i,j) = por(i,j,1)*D_o2(i,j,1)*((max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,:)) )

                     b_nh4_sub(i,j) = por(i,j,1)*D_nh4(i,j,1)*((max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,:)) )

                     b_no3_sub(i,j) = por(i,j,1)*D_no3(i,j,1)*((max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,:)) )

                     b_dic_sub(i,j) = por(i,j,1)*D_dic(i,j,1)*((max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,:)) )

                     b_odu_sub(i,j) = por(i,j,1)*D_odu(i,j,1)*((abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,:)) )

                     b_alk_sub(i,j) = por(i,j,1)*D_dic(i,j,1)*((max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) - c_talk(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) - c_talk(i,j,:)) )

                     !b_alk_sub(i,j) = sum(dz_cbed(:)*por(i,j,1:nk_cbed)*R_talk(i,j,:))  ! mol m-2 s-1 (net production of alklinity from organic matter degradation)

                     b_ca2_sub(i,j) = por(i,j,1)*D_ca2(i,j,1)*((max(0.0,fn_sw_ca2(cobalt%btm_salt(i,j), cobalt%Rho_0)) - c_ca2(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,fn_sw_ca2(cobalt%btm_salt(i,j), cobalt%Rho_0)) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,fn_sw_ca2(cobalt%btm_salt(i,j), cobalt%Rho_0)) - c_ca2(i,j,:)) )

                     b_po4_sub(i,j) = por(i,j,1)*D_po4(i,j,1)*((max(0.0,cobalt%f_po4(i,j,nk) * cobalt%Rho_0) - c_po4(i,j,1))/(dz_cbed(1)/2.0)) + &
                        por(i,j,1)*w(i,j,1)*max(0.0,cobalt%f_po4(i,j,nk) * cobalt%Rho_0) + &
                        sum(dz_cbed(:)*por(i,j,1:nk_cbed)* bioirri(i,j,:)*(max(0.0,cobalt%f_po4(i,j,nk) * cobalt%Rho_0) - c_po4(i,j,:)) )


                     ! ! check if there is any NaN or inf in b_odu. not required. can be deleted.
                     ! if (ieee_is_nan(b_odu_sub(i,j)) .or. .not. ieee_is_finite(b_odu_sub(i,j))) then
                     !    b_odu_sub(i,j) = 0.0
                     ! endif

                     cbed%cbed_b_o2_acc(i,j)  = cbed%cbed_b_o2_acc(i,j)  + b_o2_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_dic_acc(i,j)  = cbed%cbed_b_dic_acc(i,j)  + b_dic_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_nh4_acc(i,j)  = cbed%cbed_b_nh4_acc(i,j)  + b_nh4_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_no3_acc(i,j)  = cbed%cbed_b_no3_acc(i,j)  + b_no3_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_alk_acc(i,j)  = cbed%cbed_b_alk_acc(i,j)  + b_alk_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_odu_acc(i,j)  = cbed%cbed_b_odu_acc(i,j)  + b_odu_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_ca2_acc(i,j) = cbed%cbed_b_ca2_acc(i,j)  + b_ca2_sub(i,j)  * (1.0 / real(n_sub(i,j)))
                     cbed%cbed_b_po4_acc(i,j) = cbed%cbed_b_po4_acc(i,j)  + b_po4_sub(i,j)  * (1.0 / real(n_sub(i,j)))


                     !


                     ! 2. Source-sink calculations
                     ! Change `* dt` to `* dt_sub` in all these equations to update the concentrations incrementally in each sub-step, which will help prevent negative concentrations and ensure stability.

                     do k = 1, nk_cbed

                        cbed%f_o2(i,j,k)  = cbed%f_o2(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) - &
                           (2.0*R_nox(i,j,k)+R_oduox(i,j,k)) + bioirri(i,j,k)*(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,k)) )*dt_sub(i,j)

                        cbed%f_om1(i,j,k) = cbed%f_om1(i,j,k) + ( - (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k)) )*dt_sub(i,j)

                        cbed%f_om2(i,j,k) = cbed%f_om2(i,j,k) + ( - (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k)) )*dt_sub(i,j)

                        cbed%f_om3(i,j,k) = cbed%f_om3(i,j,k) + ( - (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k)) )*dt_sub(i,j)

                        cbed%f_nh4(i,j,k) = cbed%f_nh4(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                           ( - R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,k)) )*dt_sub(i,j)

                        cbed%f_no3(i,j,k) = cbed%f_no3(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
                           (R_nox(i,j,k) - 0.6*R_ana(i,j,k)) + bioirri(i,j,k)*(max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,k)) )*dt_sub(i,j)

                        cbed%f_dic(i,j,k) = cbed%f_dic(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                           1.0*svf(i,j,k)/por(i,j,k)*(R_diss_arag(i,j,k) + R_diss_calc(i,j,k)) - 1.0*(R_prec_arag(i,j,k) + R_prec_calc(i,j,k)) + &
                           bioirri(i,j,k)*(max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,k)) )*dt_sub(i,j)

                        cbed%f_odu(i,j,k) = cbed%f_odu(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
                           R_oduox(i,j,k) - svf(i,j,k)/por(i,j,k)*odu_depo(i,j,k)  + bioirri(i,j,k)*(abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,k)) )*dt_sub(i,j)

                        cbed%f_talk(i,j,k) = cbed%f_talk(i,j,k) + ( + R_talk(i,j,k) + bioirri(i,j,k)*(max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) - c_talk(i,j,k)) )*dt_sub(i,j)

                        cbed%f_calc(i,j,k) = cbed%f_calc(i,j,k) + ( - R_diss_calc(i,j,k) + por(i,j,k)/svf(i,j,k)*R_prec_calc(i,j,k) )*dt_sub(i,j)

                        cbed%f_arag(i,j,k) = cbed%f_arag(i,j,k) + ( - R_diss_arag(i,j,k) + por(i,j,k)/svf(i,j,k)*R_prec_arag(i,j,k) )*dt_sub(i,j)

                        cbed%f_ca2(i,j,k) = cbed%f_ca2(i,j,k) + ( svf(i,j,k)/por(i,j,k)*(R_diss_arag(i,j,k) + R_diss_calc(i,j,k)) - (R_prec_arag(i,j,k) + R_prec_calc(i,j,k)) + &
                           bioirri(i,j,k)*(max(0.0,fn_sw_ca2(cobalt%btm_salt(i,j), cobalt%Rho_0)) - c_ca2(i,j,k)) )*dt_sub(i,j)

                        cbed%f_po4(i,j,k) = cbed%f_po4(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*p_2_c(i,j)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                           bioirri(i,j,k)*(max(0.0,cobalt%f_po4(i,j,nk)*cobalt%Rho_0) - c_po4(i,j,k)) )*dt_sub(i,j)



                     enddo

                     ! 3. Implicit Transport
                     ! Pass `dt_sub` into vertdiff_CBED instead of `dt`
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om1, "f_om1", Db,    w, svf, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om2, "f_om2", Db,    w, svf, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om3, "f_om3", Db,    w, svf, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_o2,  "f_o2", D_o2,   w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_nh4, "f_nh4", D_nh4, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_no3, "f_no3", D_no3, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_dic, "f_dic", D_dic, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_odu, "f_odu", D_odu, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_talk, "f_talk", D_dic, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_calc, "f_calc", Db, w, svf, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_arag, "f_arag", Db, w, svf, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_ca2, "f_ca2", D_ca2, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)
                     call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_po4, "f_po4", D_po4, w, por, grid_kmt, dt_sub(i,j), tau, i,i,j,j,i,i,j,j,nk, nk_cbed)

                  enddo
                  ! --- END ADAPTIVE SUB-STEPPING LOOP ---
               endif
            enddo
         enddo



      else
         ! Source-sink calculations (no adaptive time-stepping)
         do j = jsc, jec
            do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  do k = 1, nk_cbed

                     cbed%f_o2(i,j,k)  = cbed%f_o2(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) - &
                        (2.0*R_nox(i,j,k)+R_oduox(i,j,k)) + bioirri(i,j,k)*(max(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0) - c_o2(i,j,k)) )*dt

                     cbed%f_om1(i,j,k) = cbed%f_om1(i,j,k) + ( - (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k)) )*dt

                     cbed%f_om2(i,j,k) = cbed%f_om2(i,j,k) + ( - (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k)) )*dt

                     cbed%f_om3(i,j,k) = cbed%f_om3(i,j,k) + ( - (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k)) )*dt

                     cbed%f_nh4(i,j,k) = cbed%f_nh4(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                        ( - R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(max(0.0,cobalt%f_nh4(i,j,nk)*cobalt%Rho_0) - c_nh4(i,j,k)) )*dt

                     cbed%f_no3(i,j,k) = cbed%f_no3(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
                        (R_nox(i,j,k) - 0.6*R_ana(i,j,k)) + bioirri(i,j,k)*(max(0.0,cobalt%btm_no3(i,j)*cobalt%Rho_0) - c_no3(i,j,k)) )*dt

                     cbed%f_dic(i,j,k) = cbed%f_dic(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                        bioirri(i,j,k)*(max(0.0,cobalt%btm_dic(i,j)*cobalt%Rho_0) - c_dic(i,j,k)) )*dt

                     cbed%f_odu(i,j,k) = cbed%f_odu(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
                        R_oduox(i,j,k) - svf(i,j,k)/por(i,j,k)*odu_depo(i,j,k)  + bioirri(i,j,k)*(abs(min(0.0,cobalt%btm_o2(i,j)*cobalt%Rho_0)) - c_odu(i,j,k)) )*dt

                     cbed%f_talk(i,j,k) = cbed%f_talk(i,j,k) + ( + R_talk(i,j,k) + bioirri(i,j,k)*(max(0.0,cobalt%btm_alk(i,j)*cobalt%Rho_0) - c_talk(i,j,k)) )*dt

                  enddo
               endif
            enddo
         enddo


         ! call vertdiff_CBED. This updates the fields.
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om1, "f_om1", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om2, "f_om2", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om3, "f_om3", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_o2,  "f_o2", D_o2,   w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_nh4, "f_nh4", D_nh4, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_no3, "f_no3", D_no3, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_dic, "f_dic", D_dic, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_odu, "f_odu", D_odu, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
         call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_talk, "f_talk", D_dic, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)


      endif

      ! calculate fractions of total organic matter flux assigned to each reactivity class as a func of bathymetric depth (m)
      ! writing here for testing. delete later.
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               if (cbed%use_depth_dependent_OM_frac) then
                  frac_OM1(i,j) = min(0.80, 0.65*(100.0/cobalt%zt(i,j,nk))**0.5)
                  frac_OM3(i,j) = max(0.03, 0.04*(100.0/cobalt%zt(i,j,nk))**(-0.3))
                  frac_OM2(i,j) = 1.0 - (frac_OM1(i,j) + frac_OM3(i,j))
               else
                  frac_OM1(i,j) = 0.70
                  frac_OM2(i,j) = 0.20
                  frac_OM3(i,j) = 0.10
               endif
            endif
         enddo; enddo

      ! some other diags , other local variables
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               cbed%cbed_k1(i,j) = frac_OM1(i,j)
               cbed%cbed_k2(i,j) = p_2_c(i,j)
               cbed%cbed_k3(i,j) = real(n_sub(i,j))
            endif
         enddo;enddo





      ! ! ask NIKI: if this would lead to different ans due to indexing and domain decomposition.
      ! if (cbed%do_adaptive_time_stepping) then

      !    ! --- ADAPTIVE TIME-STEPPING CALCULATION (O2, NO3, NH4, ODU) ---

      !    n_sub = 1 ! Default to 1 macro step

      !    do j = jsc, jec
      !       do i = isc, iec
      !          if (grid_kmt(i,j) > 0) then
      !             do k = 1, nk_cbed

      !                ! ---------------------------------------------------------
      !                ! 1. OXYGEN CONSTRAINT
      !                ! Sinks: Aerobic Respiration, Nitrification, ODU Oxidation
      !                ! ---------------------------------------------------------
      !                if (c_o2(i,j,k) > 1.0e-6) then
      !                   max_o2_sink = svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) + &
      !                      (2.0*R_nox(i,j,k)+R_oduox(i,j,k))

      !                   ! The maximum sink is the total amount of O2 that could be consumed in this time step based on the current concentrations and reaction rates.
      !                   ! We divide by 0.5*c_o2, to get the number of sub-steps needed to ensure that we don't drop the cell concentration more than half in one time step.
      !                   ! This is a conservative estimate to ensure we don't overshoot and get negative concentrations.
      !                   ! We then take the ceiling of this number to get the number of sub-steps needed to ensure that we don't consume more O2 than is available in any sub-step.
      !                   ! We do this for O2, NO3, and NH4 and take the maximum number of sub-steps required among the three to ensure that we don't violate any of the constraints.

      !                   n_req_o2 = ceiling( (max_o2_sink * dt) / (0.8 * c_o2(i,j,k)) )
      !                   if (n_req_o2 > n_sub) n_sub = n_req_o2
      !                endif

      !                ! ---------------------------------------------------------
      !                ! 2. NITRATE CONSTRAINT
      !                ! Sinks: Denitrification, Anammox
      !                ! ---------------------------------------------------------
      !                if (c_no3(i,j,k) > 1.0e-6) then
      !                   max_no3_sink = svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + R_ana(i,j,k)

      !                   n_req_no3 = ceiling( (max_no3_sink * dt) / (0.8 * c_no3(i,j,k)) )
      !                   if (n_req_no3 > n_sub) n_sub = n_req_no3
      !                endif

      !                ! ---------------------------------------------------------
      !                ! 3. AMMONIUM CONSTRAINT
      !                ! Sinks: Nitrification, Anammox
      !                ! ---------------------------------------------------------
      !                if (c_nh4(i,j,k) > 1.0e-6) then
      !                   max_nh4_sink = R_nox(i,j,k) + R_ana(i,j,k)

      !                   n_req_nh4 = ceiling( (max_nh4_sink * dt) / (0.8 * c_nh4(i,j,k)) )
      !                   if (n_req_nh4 > n_sub) n_sub = n_req_nh4
      !                endif

      !                ! ---------------------------------------------------------
      !                ! 4. ODU CONSTRAINT
      !                ! Sinks: ODU oxidation
      !                ! ---------------------------------------------------------
      !                if (c_odu(i,j,k) > 1.0e-6) then
      !                   max_odu_sink = R_oduox(i,j,k)

      !                   n_req_odu = ceiling( (max_odu_sink * dt) / (0.8 * c_odu(i,j,k)) )
      !                   if (n_req_odu > n_sub) n_sub = n_req_odu
      !                endif

      !             enddo
      !          endif
      !       enddo
      !    enddo

      !    ! Cap the maximum number of sub-steps to prevent the ESM from hanging
      !    ! Note: This cap can be increased depending on how aggressive the coastal fluxes may get.
      !    n_sub = min(n_sub, 60)
      !    dt_sub = dt / real(n_sub)


      !    ! --- BEGIN ADAPTIVE SUB-STEPPING LOOP ---
      !    do sub_step = 1, n_sub

      !       ! 1. Re-evaluate positive concentrations for this specific sub-step
      !       do j = jsc, jec
      !          do i = isc, iec
      !             if (grid_kmt(i,j) > 0) then
      !                do k = 1, nk_cbed
      !                   c_om1(i,j,k) = max(0.0, cbed%f_om1(i,j,k))
      !                   c_om2(i,j,k) = max(0.0, cbed%f_om2(i,j,k))
      !                   c_om3(i,j,k) = max(0.0, cbed%f_om3(i,j,k))
      !                   c_o2(i,j,k)  = max(0.0, cbed%f_o2(i,j,k))
      !                   c_no3(i,j,k) = max(0.0, cbed%f_no3(i,j,k))
      !                   c_nh4(i,j,k) = max(0.0, cbed%f_nh4(i,j,k))
      !                   c_dic(i,j,k) = max(0.0, cbed%f_dic(i,j,k))
      !                   c_odu(i,j,k) = max(0.0, cbed%f_odu(i,j,k))
      !                   c_talk(i,j,k) = max(0.0, cbed%f_talk(i,j,k))


      !                   ! O₂ reaction rates
      !                   R_om1_o2(i,j,k) = k1(i,j)*c_om1(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   R_om2_o2(i,j,k) = k2(i,j)*c_om2(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   R_om3_o2(i,j,k) = k3(i,j)*c_om3(i,j,k)*(c_o2(i,j,k)/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   ! NO₃ reaction rates
      !                   R_om1_no3(i,j,k) = k_adj_denit*k1(i,j)*c_om1(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   R_om2_no3(i,j,k) = k_adj_denit*k2(i,j)*c_om2(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   R_om3_no3(i,j,k) = k_adj_denit*k3(i,j)*c_om3(i,j,k)*(c_no3(i,j,k)/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   ! ODU reaction rates
      !                   R_om1_anoxic(i,j,k) = k_adj_anoxia*k1(i,j)*c_om1(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   R_om2_anoxic(i,j,k) = k_adj_anoxia*k2(i,j)*c_om2(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)
      !                   R_om3_anoxic(i,j,k) = k_adj_anoxia*k3(i,j)*c_om3(i,j,k)*(ks_no3/(ks_no3 + c_no3(i,j,k)))*(ks_o2/(ks_o2 + c_o2(i,j,k))) * Q10_factor(i,j)

      !                   ! dic
      !                   R_dic_om1(i,j,k) = (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k))
      !                   R_dic_om2(i,j,k) = (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k))
      !                   R_dic_om3(i,j,k) = (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k))
      !                   ! nitrification
      !                   R_nox(i,j,k) = k_nox*c_nh4(i,j,k)*c_o2(i,j,k) * Q10_factor(i,j)
      !                   ! anammox
      !                   R_ana(i,j,k) = k_ana*c_nh4(i,j,k)*c_no3(i,j,k) * Q10_factor(i,j) !* (ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
      !                   ! ODU oxidation
      !                   R_oduox(i,j,k) = k_oduox*c_odu(i,j,k)*c_o2(i,j,k) * Q10_factor(i,j)
      !                   odu_depo(i,j,k) = (R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))*min(1.0, 0.233*(w(i,j,k)*100.0*spery)**0.336)

      !                   ! TA calculation
      !                   R_talk(i,j,k) = svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) + &
      !                      svf(i,j,k)/por(i,j,k)*(0.8+1.0/cobalt%c_2_n)*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
      !                      svf(i,j,k)/por(i,j,k)*(1.0+1.0/cobalt%c_2_n)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
      !                      2.0*R_nox(i,j,k) - 1.0*R_oduox(i,j,k)

      !                enddo
      !             endif
      !          enddo
      !       enddo


      !       ! 2. Source-sink calculations
      !       ! Change `* dt` to `* dt_sub` in all these equations to update the concentrations incrementally in each sub-step, which will help prevent negative concentrations and ensure stability.
      !       do j = jsc, jec
      !          do i = isc, iec
      !             if (grid_kmt(i,j) > 0) then
      !                do k = 1, nk_cbed

      !                   cbed%f_o2(i,j,k)  = cbed%f_o2(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) - &
      !                      (2.0*R_nox(i,j,k)+R_oduox(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_o2(i,j)*cobalt%Rho_0 - c_o2(i,j,k)) )*dt_sub

      !                   cbed%f_om1(i,j,k) = cbed%f_om1(i,j,k) + ( - (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k)) )*dt_sub

      !                   cbed%f_om2(i,j,k) = cbed%f_om2(i,j,k) + ( - (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k)) )*dt_sub

      !                   cbed%f_om3(i,j,k) = cbed%f_om3(i,j,k) + ( - (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k)) )*dt_sub

      !                   cbed%f_nh4(i,j,k) = cbed%f_nh4(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
      !                      ( - R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%f_nh4(i,j,nk)*cobalt%Rho_0 - c_nh4(i,j,k)) )*dt_sub

      !                   cbed%f_no3(i,j,k) = cbed%f_no3(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
      !                      (R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_no3(i,j)*cobalt%Rho_0 - c_no3(i,j,k)) )*dt_sub

      !                   cbed%f_dic(i,j,k) = cbed%f_dic(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
      !                      bioirri(i,j,k)*(cobalt%btm_dic(i,j)*cobalt%Rho_0 - c_dic(i,j,k)) )*dt_sub

      !                   cbed%f_odu(i,j,k) = cbed%f_odu(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
      !                      R_oduox(i,j,k) - svf(i,j,k)/por(i,j,k)*odu_depo(i,j,k)  + bioirri(i,j,k)*(0.0 - c_odu(i,j,k)) )*dt_sub

      !                   cbed%f_talk(i,j,k) = cbed%f_talk(i,j,k) + ( + R_talk(i,j,k) + bioirri(i,j,k)*(cobalt%btm_alk(i,j)*cobalt%Rho_0 - c_talk(i,j,k)) )*dt_sub

      !                enddo
      !             endif
      !          enddo
      !       enddo

      !       ! 3. Implicit Transport
      !       ! Pass `dt_sub` into vertdiff_CBED instead of `dt`
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om1, "f_om1", Db,    w, svf, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om2, "f_om2", Db,    w, svf, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om3, "f_om3", Db,    w, svf, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_o2,  "f_o2", D_o2,   w, por, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_nh4, "f_nh4", D_nh4, w, por, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_no3, "f_no3", D_no3, w, por, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_dic, "f_dic", D_dic, w, por, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_odu, "f_odu", D_odu, w, por, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !       call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_talk, "f_talk", D_dic, w, por, grid_kmt, dt_sub, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)

      !    enddo
      !    ! --- END ADAPTIVE SUB-STEPPING LOOP ---



      ! else
      !    ! Source-sink calculations
      !    do j = jsc, jec; do i = isc, iec  !{
      !          if (grid_kmt(i,j) .gt. 0) then
      !             do k = 1, nk_cbed

      !                !cbed%f_tr1(i,j,k) = cbed%f_tr1(i,j,k) + 0.01 * k !fictitious dubious dynamics for testing purposes

      !                cbed%f_o2(i,j,k)  = cbed%f_o2(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) - &
      !                   (2.0*R_nox(i,j,k)+R_oduox(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_o2(i,j)*cobalt%Rho_0 - c_o2(i,j,k)) )*dt

      !                cbed%f_om1(i,j,k) = cbed%f_om1(i,j,k) + ( - (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k)) )*dt

      !                cbed%f_om2(i,j,k) = cbed%f_om2(i,j,k) + ( - (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k)) )*dt

      !                cbed%f_om3(i,j,k) = cbed%f_om3(i,j,k) + ( - (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k)) )*dt

      !                cbed%f_nh4(i,j,k) = cbed%f_nh4(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
      !                   ( - R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%f_nh4(i,j,nk)*cobalt%Rho_0 - c_nh4(i,j,k)) )*dt

      !                cbed%f_no3(i,j,k) = cbed%f_no3(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
      !                   (R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_no3(i,j)*cobalt%Rho_0 - c_no3(i,j,k)) )*dt

      !                cbed%f_dic(i,j,k) = cbed%f_dic(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
      !                   bioirri(i,j,k)*(cobalt%btm_dic(i,j)*cobalt%Rho_0 - c_dic(i,j,k)) )*dt

      !                cbed%f_odu(i,j,k) = cbed%f_odu(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
      !                   R_oduox(i,j,k) - svf(i,j,k)/por(i,j,k)*odu_depo(i,j,k)  + bioirri(i,j,k)*(0.0 - c_odu(i,j,k)) )*dt

      !                cbed%f_talk(i,j,k) = cbed%f_talk(i,j,k) + ( + R_talk(i,j,k) + bioirri(i,j,k)*(cobalt%btm_alk(i,j)*cobalt%Rho_0 - c_talk(i,j,k)) )*dt


      !             enddo
      !          endif
      !       enddo;enddo


      !    ! call vertdiff_CBED. This updates the fields.
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om1, "f_om1", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om2, "f_om2", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om3, "f_om3", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_o2,  "f_o2", D_o2,   w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_nh4, "f_nh4", D_nh4, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_no3, "f_no3", D_no3, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_dic, "f_dic", D_dic, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_odu, "f_odu", D_odu, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      !    call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_talk, "f_talk", D_dic, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)


      ! endif

      !------------------
      ! Write checksum
      !------------------
      stdoutunit=stdout()
      !write(stdoutunit,'("cbedsums: ",Z16.16,ES24.16)') sum(cbed%f_o2(:,:,:)),sum(cbed%f_dic(:,:,:)),sum(cbed%f_om2(:,:,:))

      ! Print Oxygen (f_o2) sum
      write(stdoutunit,'("cbedsum f_o2:  ",Z16.16,1X,ES24.16)') &
         sum(cbed%f_o2(:,:,:)), sum(cbed%f_o2(:,:,:))

      ! Print Dissolved Inorganic Carbon (f_dic) sum
      write(stdoutunit,'("cbedsum f_dic: ",Z16.16,1X,ES24.16)') &
         sum(cbed%f_dic(:,:,:)), sum(cbed%f_dic(:,:,:))

      ! Print Organic Matter (f_om2) sum
      write(stdoutunit,'("cbedsum f_om2: ",Z16.16,1X,ES24.16)') &
         sum(cbed%f_om2(:,:,:)), sum(cbed%f_om2(:,:,:))




      !!==================================================================================================================
      !!The rest of this subrouine that follows is a copy of the COBALT code.
      !!It must be replaced by CBED calculations for 'btf' fluxes which are "set" for COBALT at the end of this subroutine.
      !!==================================================================================================================

      ! Calculate the bottom conditions and the fluxes to the bottom for diagnostics and benthic flux calculations.
      ! MOM4/5 used the bottom grid cell, but MOM6 often has a number of vanishingly thin layers overlying the bottom.
      ! Grid scale noise in these layers can occur, particularly for quantities with large bottom fluxes.  COBALT thus
      ! uses conditions over a specified bottom layer thickness (cobalt%bottom_thickness, default = 1m) for bottom calcs.

      do j = jsc, jec; do i = isc, iec  !{
            if (grid_kmt(i,j) .gt. 0) then !{


               ! Calculate the processing of organic matter in the sediment.  The fate of organic matter is partitioned
               ! between burial (i.e., removal from the system), aerobic remineralization, remineralization via
               ! denitrification, and remineralization via sulfate reduction.  Note that the latter pathway is effectively
               ! a "catch all" for any other anaerobic pathway and the sulfate cycle is not explicitly modeled.
               k = grid_kmt(i,j)
               if (cobalt%fntot_btm(i,j) .gt. 0.0) then !{

                  ! The Burial flux estimates are based on Dunne et al., 2007. A synthesis of global particle export from
                  ! the surface ocean and cycling through the ocean interior and on the seafloor.  Global Biogeochemical
                  ! Cycles. Vol. 21, GB4006, doi:10.1029/2006GB002907.  See Figure 2, eq. (3).  The default units of this
                  ! relationship are mmoles C m-2 day-1, and the local variable "fpoc_btm" is used to create a bottom flux
                  ! in these units.
                  !
                  ! As described in Dunne et al., (2007) this relationship was generally developed for deeper ocean areas
                  ! and its validity in shallow areas is unclear.  Past experiments suggest that it may overestimate burial
                  ! in shallow areas, resulting in large nutrient losses that are inconsistent with observations.  The
                  ! parameter "z_burial" thus provides a depth scale (an effective "half-saturation") for ramping up burial
                  ! from 0 to its full value.
                  !
                  ! Since burial is highly uncertain and often used in global earth system simulations to balance inputs and
                  ! outputs, a dimensionless scaling factor (cobalt%scale_burial) has also been included.
                  fpoc_btm = cobalt%fntot_btm(i,j)*cobalt%c_2_n*sperd*1000.0
                  cobalt%frac_burial(i,j) = 0.013 + 0.53*fpoc_btm**2.0/((7.0+fpoc_btm)**2.0) * &
                     cobalt%zt(i,j,k) / (cobalt%z_burial + cobalt%zt(i,j,k))
                  cobalt%frac_burial(i,j) = cobalt%scale_burial*cobalt%frac_burial(i,j)
                  cobalt%fn_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%fntot_btm(i,j)
                  cobalt%fp_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%fptot_btm(i,j)

                  !-----------------------------------------------
                  !!! replace cobalt burial with CBED burial calculation
                  !-----------------------------------------------

                  !cobalt%frac_burial(i,j) = cbed_burial_frac(i,j)
                  cobalt%fn_burial(i,j) = cbed%burial_om(i,j) * (1.0/cobalt%c_2_n)
                  cobalt%fp_burial(i,j) = min(1.0, cbed_burial_frac(i,j)) * cobalt%fptot_btm(i,j)

                  !!!------------------------------------------

                  ! Denitrification follows Middelburg et al., 1996. Denitrification in marine sediments: a modeling study
                  ! Global Biogeochemical Cycles 10(4).  pp. 661-673.  https://doi.org/10.1029/96GB02562. COBALT uses the
                  ! carbon flux-based relationship based on Middelburg's first extraction of his metamodel (the first
                  ! equation in Section 3.4 of the paper).  This relationship requires a flux to the benthos in micromoles C
                  ! cm-2 day-1.  This means that fpoc_btm defined for the burial calculation above must be multiplied by:
                  !
                  ! 1e3 micromoles/millimole*1e-4 cm2/m2 = 0.1
                  !
                  ! to get the proper units.  The Middelburg relationship yields a rate at which arriving particulate organic
                  ! carbon is denitrified in micromoles C cm-2 day-1.  This is converted to a rate at which arriving
                  ! particulate organic nitrogen denitrified in moles N m-2 sec-1 by dividing by:
                  !
                  ! c_2_n*sperd*1e6 micromoles/mole*1e-4 cm2/m2 = c_2_n*sperd*100
                  !
                  ! The nitrate demand associated with this denitrification (fno3denit_sed) is obtained by multiplying the
                  ! resulting value by the moles of NO3 required to denitrify each mole of organic N (n_2_n_denit).
                  !
                  ! A number of limiters are applied to support global application.  First, the C flux used in the
                  ! Middelburg relationship is capped at 43.0 micromoles C cm-2 day-1 to avoid anomalous extrapolation.
                  ! Second, denitrification is slowed when bottom nitrate is low by a) scaling rates with a nitrate
                  ! half-saturation constant with (k_no3_denit), b) preventing the exhaustion of bottom nitrate over
                  ! single time step, and c) limiting the total amount of organic carbon denitrified to that arriving at
                  ! the sediment minus that which was buried. Finally, to prevent excessive denitrification in very shallow
                  ! areas, a depth scale (z_denit) was included to ramp up rates to full Middelburg values only in deeper
                  ! waters.
                  log10_fpoc_btm = log10(min(43.0,0.1*fpoc_btm))
                  cobalt%fno3denit_sed(i,j) = min(cobalt%btm_no3(i,j)*cobalt%bottom_thickness*cobalt%Rho_0*r_dt,  &
                     min((cobalt%fntot_btm(i,j)-cobalt%fn_burial(i,j))*cobalt%n_2_n_denit, &
                     10.0**(-0.9543+0.7662*log10_fpoc_btm - 0.235*log10_fpoc_btm**2.0)/(cobalt%c_2_n*sperd*100.0)* &
                     cobalt%n_2_n_denit*cobalt%btm_no3(i,j)/(cobalt%k_no3_denit + cobalt%btm_no3(i,j)))) * &
                     cobalt%zt(i,j,k) / (cobalt%z_denit + cobalt%zt(i,j,k))

                  ! Calculate the rate of organic matter degradation in the sediment after accounting for burial
                  ! and denitrification.  Two pathways are tracked:
                  !
                  ! fnoxic_sed (moles N m-2 sec-1) accounts for organic material remineralized by processes using oxygen
                  ! *within the sediment*, resulting in an oxygen demand at the sediment-water interface.  These
                  ! include direct aerobic remineralization and sulfate reduction/HS- oxidation (see stoichiometry
                  ! for details).  Note that the partitioning between these two pathways is not calculated, just the
                  ! combined effect.
                  !
                  ! fnso4_sed (moles N m-2 sec-1) accounts for organic material that only undergoes sulfate reduction in
                  ! the sediment, but not HS- oxidation.  This results in HS- released from the sediment.  The latent O2
                  ! demand from the HS- is tracked if "do_fnso4red_sed = true".  The amount of material falling into this
                  ! category is equal to the remainder after all other pathways are accounted for and it generally only
                  ! occurs under anoxic conditions in COBALT.  The added O2 demand from HS- can result in negative O2
                  ! concentrations that should be interpreted as 0 moles O2 kg-1 + additional O2 demand from HS-
                  !
                  ! NOTE: The O2 demand from the sediment is calculated by multiplying the organic matter remineralized
                  !       by the O2 per N (i.e., fnoxic_sed*o2_2_nh4 or (fnoxic_sed+fnso4red_sed)*o2_2_nh4 if
                  !       do_fnso4red_sed is true).
                  ! NOTE: fnso4red_sed is not the total sulfate reduction in the sediment, only that which is not paired
                  !       with subsequent HS- oxidation.  COBALT does not calculate the total sulfate reduction in seds.
                  ! NOTE: The maximum organic remin supported by local O2 is:
                  ! btm_o2(moles O2 kg-1)*bottom_thickness(m)*density(kg m-3)* 1/dt(s-1)*molN/molO2 = moles N m-2 s-1
                  !
                  ! The thickness of the bottom boundary layer (cobalt%bottom_thickness) impacts this upper bound.
                  ! Efforts are underway to implement a more dynamic bottom boundary layer scheme.
                  !
                  if (cobalt%btm_o2(i,j) .gt. cobalt%o2_min) then  !{
                     cobalt%fnoxic_sed(i,j) = max(0.0, min(cobalt%btm_o2(i,j)*cobalt%bottom_thickness* &
                        cobalt%Rho_0*r_dt*(1.0/cobalt%o2_2_nh4), &
                        cobalt%fntot_btm(i,j) - cobalt%fn_burial(i,j) - &
                        cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit))
                  else
                     cobalt%fnoxic_sed(i,j) = 0.0
                  endif !}
                  cobalt%fnso4red_sed(i,j) = max(0.0, cobalt%fntot_btm(i,j)-cobalt%fnoxic_sed(i,j)- &
                     cobalt%fn_burial(i,j)-cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)
               else
                  cobalt%fnso4red_sed(i,j) = 0.0
                  cobalt%fno3denit_sed(i,j) = 0.0
                  cobalt%fnoxic_sed(i,j) = 0.0
               endif !}

               !
               ! Iron flux from the sediment
               !

               ! Iron from sediment (Dale, 2015).  The maximum release from the sediment is set by ffe_sed_max.  The
               ! hyperbolic tangent requires the flux of carbon to the sediments (as mmoles m-2 day-1) in the numerator
               ! and the bottom water oxygen concentration (in microMolar units) in the denominator. Note that ffe_sed_max
               ! was converted to moles Fe m-2 sec-1 during parameter input, so ffe_sed is in moles Fe m-2 sec-1
               cobalt%ffe_sed(i,j) = cobalt%ffe_sed_max * tanh( (cobalt%fntot_btm(i,j)*cobalt%c_2_n*sperd*1.0e3)/ &
                  max(cobalt%btm_o2(i,j)*1.0e6,epsln) )

               ! Additional coastal iron (Optional, default fe_coast = 0)
               !
               ! Coarse resolution models and/or intermediate resolution models in areas with exceptionally steep bathymetry
               ! can under-represent coastal iron because they don't resolve shallow regions. An option to add iron through
               ! the vertical face of the land mass has thus been included.  The flux is posed as a fraction (fe_coast) of
               ! the sediment Fe flux (moles Fe m-2 sec-1) that would have resulted from the sinking organic matter flux and
               ! O2 level of the adjacent waters.  This is then spread across the layer mass (rho_dzt(i,j,k)) to give an input
               ! in moles Fe kg-1 sec-1. Conceptually, this can be thought of as a net iron flux resulting from the fraction
               ! of the sinking flux that would have been intercepted at shallower depths were the model resolution finer.
               ! The default value of fe_coast is 0 (i.e., only the explicitly resolved benthic flux is included).
               !
               ! Old Expression:
               ! cobalt%jfe_coast(i,j,1) = cobalt%fe_coast * mask_coast(i,j) * grid_tmask(i,j,1) / &
               !     sqrt(grid_dat(i,j))
               !
               do k = 1, nk !{
                  if (cobalt%fe_coast == 0.0) then
                     cobalt%jfe_coast(i,j,k) = 0.0
                  else
                     cobalt%jfe_coast(i,j,k) = cobalt%fe_coast*dzt(i,j,k)*mask_coast(i,j)*grid_tmask(i,j,k)* &
                        cobalt%ffe_sed_max*tanh( ( (cobalt%f_ndet(i,j,k)*cobalt%wsink+ &
                        phyto(SMALL)%f_n(i,j,k)*phyto(SMALL)%vmove(i,j,k)+ &
                        phyto(MEDIUM)%f_n(i,j,k)*phyto(MEDIUM)%vmove(i,j,k)+ &
                        phyto(LARGE)%f_n(i,j,k)*phyto(LARGE)%vmove(i,j,k)+ &
                        phyto(DIAZO)%f_n(i,j,k)*phyto(DIAZO)%vmove(i,j,k))*cobalt%c_2_n*sperd*1.0e3 )/ &
                        max(cobalt%f_o2(i,j,k)*1.0e6,epsln) )/rho_dzt(i,j,k)
                  endif
               enddo  !} k

               ! Have ffe_geotherm default to zero if the internal_heat variable
               ! needed to calculate it is not available (if geothermal heating is disabled).
               if(present(internal_heat)) then
                  cobalt%ffe_geotherm(i,j) = cobalt%ffe_geotherm_ratio*internal_heat(i,j)*4184.0/dt
               else
                  cobalt%ffe_geotherm(i,j) = 0.0
               endif

               !
               ! Calcium carbonate flux and burial, based on Dunne et al., 2012
               !
               ! phi_surfresp_cased = 0.14307   ! const for enhanced diss., surf sed respiration (dimensionless)
               ! phi_deepresp_cased = 4.1228    ! const for enhanced diss., deep sed respiration (dimensionless)
               ! alpha_cased = 2.7488 ! exponent controlling non-linearity of deep dissolution
               ! beta_cased = -2.2185 ! exponent controlling non-linearity of effective thickness
               ! gamma_cased = 0.03607/spery   ! dissolution rate constant
               ! Co_cased = 8.1e3        ! moles CaCo3 m-3 for pure calcite sediment with porosity = 0.7
               !
               ! if cased_steady is true, burial is calculated from Dunne's eq. (2) assuming dcased/dt = 0.
               ! This ensures that all the calcite bottom flux is partitioned between burial and redissolution.
               ! The steady state cased value of cased is calculated to reflect the changing bottom conditions.
               ! This influences the the partitioning of burial and redissolution over time, but there are
               ! no alkalinity changes/drifts associated with the long-term evolution of cased
               !
               ! If cased_steady is false, calcite is partitioned between dissolution, burial and evolving
               ! cased as described in Dunne et al. (2012).  The multi-century scale evolution of cased
               ! impacts alkalinity, but care must to ensure that cased starts in equilibrium with the
               ! mean ocean state to avoid unrealistic drifts.

               k = grid_kmt(i,j)

               ! Enhanced dissolution by fast respiration near the sediment surface, proportional
               ! to organic flux, moles Ca m-2 s-1, limited to a max 1/2 the instantaneous calcite flux
               cobalt%fcased_redis_surfresp(i,j)=min(0.5*cobalt%f_cadet_calc_btf(i,j,1), &
                  cobalt%phi_surfresp_cased*cobalt%fntot_btm(i,j)*cobalt%c_2_n)

               ! Ca-specific dissolution coeficient, depends on calcite saturation state and is enhanced by
               ! respiration deep in the sediment (s-1), non-linearity controlled by alpha_cased
               cobalt%cased_redis_coef(i,j) = cobalt%gamma_cased*max(0.0,1.0-cobalt%btm_omega_calc(i,j)+ &
                  cobalt%phi_deepresp_cased*cobalt%fntot_btm(i,j)*cobalt%c_2_n*spery)**cobalt%alpha_cased

               ! Effective thickness term that enhances burial of calcite when total sediment accumulation is high
               ! dimensionless value between 0 and 1
               cobalt%cased_redis_delz(i,j) = max(1.0, &
                  cobalt%f_lithdet_btf(i,j,1)*spery+cobalt%f_cadet_calc_btf(i,j,1)*100.0*spery)**cobalt%beta_cased

               ! calculate the sediment redissolution rate (moles Ca m-2 sec-1). This calculation is subject to
               ! three limiters: a) a maximum of 1/2 of the total cased over one time step; b) a maximum of 0.01
               ! moles Ca per day; and c) a minimum of 0.0
               cobalt%fcased_redis(i,j) = max(0.0, min(0.01/sperd, min(0.5*cobalt%f_cased(i,j,1)*r_dt,  &
                  cobalt%fcased_redis_surfresp(i,j)+cobalt%cased_redis_coef(i,j)*cobalt%cased_redis_delz(i,j)*cobalt%f_cased(i,j,1))) )

               !
               ! Old expression
               !
               !cobalt%fcased_redis(i,j) = max(0.0, min(0.01/sperd,min(0.5 * cobalt%f_cased(i,j,1) * r_dt, min(0.5 *       &
               !   cobalt%f_cadet_calc_btf(i,j,1), 0.14307 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n) +        &
               !   0.03607 / spery * max(0.0, 1.0 - cobalt%omega_calc(i,j,k) +   &
               !   4.1228 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n * spery)**(2.7488) *                        &
               !   max(1.0, cobalt%f_lithdet_btf(i,j,1) * spery + cobalt%f_cadet_calc_btf(i,j,1) * 100.0 *  &
               !   spery)**(-2.2185) * cobalt%f_cased(i,j,1))))*grid_tmask(i,j,k)

               if (cobalt%cased_steady) then
                  cobalt%fcased_burial(i,j) = cobalt%f_cadet_calc_btf(i,j,1) - cobalt%fcased_redis(i,j)
                  cobalt%f_cased(i,j,1) = cobalt%fcased_burial(i,j)*cobalt%Co_cased/cobalt%f_cadet_calc_btf(i,j,1)
               else
                  cobalt%fcased_burial(i,j) = max(0.0, cobalt%f_cadet_calc_btf(i,j,1) * cobalt%f_cased(i,j,1) / &
                     cobalt%Co_cased)
                  cobalt%f_cased(i,j,1) = cobalt%f_cased(i,j,1) + (cobalt%f_cadet_calc_btf(i,j,1) -            &
                     cobalt%fcased_redis(i,j) - cobalt%fcased_burial(i,j)) / cobalt%z_sed * dt *                &
                     grid_tmask(i,j,k)
               endif

               !
               ! Bottom flux boundaries passed to the vertical mixing routine
               ! (negative values are fluxes into the ocean)
               !
               cobalt%b_dic(i,j) =  - cobalt%fcased_redis(i,j) - cobalt%f_cadet_arag_btf(i,j,1) -       &
                  (cobalt%fntot_btm(i,j) - cobalt%fn_burial(i,j)) * cobalt%c_2_n
               cobalt%b_fed(i,j) = - cobalt%ffe_sed(i,j) - cobalt%ffe_geotherm(i,j)
               cobalt%b_nh4(i,j) = - cobalt%fntot_btm(i,j) + cobalt%fn_burial(i,j)
               cobalt%b_no3(i,j) = cobalt%fno3denit_sed(i,j)
               ! Include latent O2 demand and alkalinity effects of HS- (see stoichiometry)
               if (cobalt%do_fnso4red_sed) then
                  cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * (cobalt%fnoxic_sed(i,j) + cobalt%fnso4red_sed(i,j))
                  cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
                     cobalt%fnoxic_sed(i,j) - cobalt%fno3denit_sed(i,j)*cobalt%alk_2_n_denit - cobalt%fnso4red_sed(i,j)
               else
                  cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * cobalt%fnoxic_sed(i,j)
                  cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
                     cobalt%fnoxic_sed(i,j) - cobalt%fno3denit_sed(i,j)*cobalt%alk_2_n_denit
               endif
               cobalt%b_po4(i,j) = - cobalt%fptot_btm(i,j) + cobalt%fp_burial(i,j)
               cobalt%b_sio4(i,j)= - cobalt%fsitot_btm(i,j)

            endif !}
         enddo; enddo  !} i, j

      do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
               cobalt%f_cased(i,j,k) = 0.0
            enddo; enddo ; enddo  !} i,j,k



      ! ! set the cobalt%b_* terms. These are used in some other places in COBALT as well. So just "b_*" might not work.
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then

               ! the b_dic and b_alk are defined here such that it takes organic part from CBED and CaCO3 part from COBALT.
               ! For b_dic the diffusive flux is added as the organic part as determined by CBED porewater DIC.
               ! Later, once CaCO3 is implemented in CBED, the b_dic and b_alk can be fully calculated from CBED from the diffusive
               ! gradients of porewater DIC and alk, and COBALT terms can be removed.
               ! for the b_alk, the net production of alkalinity from organics is added to the COBALT as porewater alkalinity
               ! profiles and resulting diffive fluxes are imcomplete without CaCO3 implemented in CBED.

               cobalt%b_dic(i,j) = b_dic(i,j)
               !cobalt%b_dic(i,j) =  - cobalt%fcased_redis(i,j) - cobalt%f_cadet_arag_btf(i,j,1) +       &
               !   b_dic(i,j)

               cobalt%b_alk(i,j) = b_alk(i,j)
               !cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
               !   b_alk(i,j) !cbed_org_alk(i,j)

               ! Add the ODU flux as added oxygen demand by the sediment because released ODU will be consummed in the bottom water.
               ! In absence of BW O2, it will create -ve O2 conc in BW. cbed%odu_flux(i,j) value is negative meaning efflux of ODU from sediment.
               ! Multiply with - sign will convert it to +ve meaning it will effectively "increase" b_o2 i.e. benthic oxygen demand.

               cobalt%b_o2(i,j)  = b_o2(i,j) + (- b_odu(i,j))


               cobalt%b_nh4(i,j) = b_nh4(i,j)
               cobalt%b_no3(i,j) = b_no3(i,j)

               ! other b terms are also calculated from CBED information. These are done by replacing cobalt%burial_frac with cbed_burial_frac where they are defined.


            endif
         enddo; enddo

      ! CBED is considered in cobalt% b_alk, b_dic, b_nh4, b_no3, b_o2, and b_po4.

      call g_tracer_set_values(cobalt_tracer_list,'alk',  'btf', cobalt%b_alk ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'dic',  'btf', cobalt%b_dic ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'fed',  'btf', cobalt%b_fed ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'nh4',  'btf', cobalt%b_nh4 ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'no3',  'btf', cobalt%b_no3 ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'o2',   'btf', cobalt%b_o2  ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'po4',  'btf', cobalt%b_po4 ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'sio4', 'btf', cobalt%b_sio4,isd,jsd)

   end subroutine generic_CBED_update_from_source

end module generic_CBED
