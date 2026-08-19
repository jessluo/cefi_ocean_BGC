!> Opens the CBED runtime parameter file(s) and returns a handle that
!! generic_CBED can hand to the MOM6 get_param interface.
module cbed_param_doc

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : param_file_type, open_param_file
use MOM_io,            only : file_exists, close_file, slasher, ensembler
use MOM_io,            only : open_namelist_file, check_nml_error
!use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time_type
!use MOM_time_manager,  only : operator(+), operator(-), operator(>)

implicit none ; private

public get_CBED_param_file

!> Container for paths and parameter file names.
!type, public :: directories
!  character(len=240) :: &
!    restart_input_dir = ' ',& !< The directory to read restart and input files.
!    restart_output_dir = ' ',&!< The directory into which to write restart files.
!    output_directory = ' '    !< The directory to use to write the model output.
!  character(len=2048) :: &
!    input_filename  = ' '     !< A string that indicates the input files or how
!                              !! the run segment should be started.
!end type directories

contains
!! This subroutine is analagous to get_COBALT_param_file in the cobalt_param_doc module,
!> which is in turn analagous to get_mom_input within MOM_get_input. The additional
!> arguments used in the MOM version are commented out rather than deleted, so that they
!> are easy to restore as the CBED parameter files grow in complexity.
!>
!> The parameter file(s) are named by the parameter_filename entry of the &cbed_input_nml
!> namelist group in input.nml, e.g.
!>
!>   &cbed_input_nml
!>       parameter_filename = 'CBED_input'
!>   /
!>
!> If that group is absent, 'CBED_input' is assumed. An empty CBED_input is valid and
!> leaves every CBED parameter at its default.
subroutine get_CBED_param_file(param_file)
  type(param_file_type), optional, intent(out) :: param_file   !< A structure to parse for run-time parameters.
  !type(directories),     optional, intent(out) :: dirs         !< Container for paths and parameter filenames.
  !logical,               optional, intent(in)  :: check_params !< If present and False will stop error checking for
  !                                                             !! run-time parameters.
  !character(len=*),      optional, intent(in)  :: default_input_filename !< If present, is the value assumed for
  !                                                             !! input_filename if input_filename is no listed
  !                                                             !! in the namelist CBED_input_nml.
  ! Local variables
  integer, parameter :: npf = 5 ! Maximum number of parameter files

  character(len=240) :: &
    output_directory = ' ', &      ! Directory to use to write the model output.
    parameter_filename(npf) = ' '  ! List of files containing parameters.

  character(len=2048) :: &
    input_filename             ! A string that indicates the input files or how
                               ! the run segment should be started.
  character(len=240) :: output_dir
  integer :: unit, io, ierr, valid_param_files

  namelist /cbed_input_nml/ parameter_filename

  ! Open namelist
  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call MOM_error(FATAL,'Required namelist file input.nml does not exist.')
  endif

  ! Assume 'CBED_input' when &cbed_input_nml is absent from input.nml. This is assigned
  ! before the read so that a namelist that names only some of the npf slots still leaves
  ! the default in slot 1.
  parameter_filename(:) = ' '
  parameter_filename(1) = 'CBED_input'

  ! Read namelist parameters
  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=cbed_input_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'cbed_input_nml')
  enddo
10 call close_file(unit)

  output_dir = trim(slasher(ensembler(output_directory)))
  valid_param_files = 0
  do io=1,npf
     if (len_trim(trim(parameter_filename(io))) > 0) then
       ! open_param_file issues its own FATAL for a missing file, but the message there
       ! gives no hint about CBED or about which namelist named the file, so check here.
       if (.not. file_exists(trim(parameter_filename(io)))) then
         call MOM_error(FATAL, "get_CBED_param_file: the CBED parameter file '"// &
             trim(parameter_filename(io))//"' does not exist. When do_CBED=.true. in "// &
             "&generic_COBALT_nml, every file named by parameter_filename in "// &
             "&cbed_input_nml in input.nml must exist. An empty file is valid and "// &
             "leaves all CBED parameters at their default values.")
       endif
       call open_param_file(trim(parameter_filename(io)), param_file, &
                            component="CBED", doc_file_dir=output_dir)
        valid_param_files = valid_param_files + 1
     endif
   enddo
   if (valid_param_files == 0) call MOM_error(FATAL, "There must be at "//&
      "least 1 valid entry in parameter_filename in cbed_input_nml in input.nml.")

end subroutine get_CBED_param_file

end module cbed_param_doc
