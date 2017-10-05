#
# Second-level makefile for ONETEP
#
# This version for multiple architectures by Peter Haynes
#
# System-specific changes should be made to the appropriate
#  conf.arch files
#

EXECDIR = $(ROOTDIR)/bin
INCPATH = -I$(ROOTDIR)/include
LIBPATH = -L$(ROOTDIR)/lib/$(ARCH)

include $(ROOTDIR)/config/conf.$(ARCH)

EXTFLAGS = $(OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .o .F90

#CW
# optics_mod.o was moved ahead of hubbard_build_mod.o for reasons
# related to the DMFT interface
# eigenstates_mod.o added to ngwf_cg 
# kernel_diis_mod.F90 added to hubbard_build_mod.o
# eigenstates_mod.o added to hamiltonian_mod.o
# added hubbard_build_mod.o hubbard_init_mod.o to bandstructure
# added density_mod.o to hubbard_build_mod.o  
#END CW

#
# List of source files
#

SOURCES = atom_mod.F90 augmentation_mod.F90 bandstructure_mod.F90 \
          basis_mod.F90 cdft_cg_mod.F90 cell_mod.F90 \
          cell_grid_mod.F90 chebyshev_rep_mod.F90 classical_pot_mod.F90 \
          comms_mod.F90 \
          conduction_mod.F90 conduction_properties_mod.F90 \
          constants_mod.F90 cutoff_coulomb_mod.F90 \
          dense_mod.F90 density_mod.F90 electronic_mod.F90 eigenstates_mod.F90 \
          energy_and_force_mod.F90 enthalpy_mod.F90 esdf_key_mod.F90 \
          esdf_mod.F90 ewald_mod.F90 finite_differences_mod.F90 \
          forces_mod.F90 fourier_mod.F90 function_basis_mod.F90 \
          gaunt_coeff_mod.F90 geometry_deloc_algor.F90 \
          geometry_deloc_utils.F90 geometry_mod.F90 \
          geometry_optimiser_mod.F90 hamiltonian_mod.F90 \
          hartree_mod.F90 hf_exchange_mod.F90 optics_mod.F90 hubbard_build_mod.F90 \
          hubbard_init_mod.F90 integrals_mod.F90 ion_mod.F90 \
          is_poisson_mod.F90 is_smeared_ions_mod.F90 \
          is_solvation_mod.F90 kernel_diis_mod.F90 kernel_mod.F90 \
          kinetic_mod.F90 lnv_mod.F90 md_mod.F90 md_ions_mod.F90 \
          md_thermostat_mod.F90 multigrid_mod.F90 multigrid_methods_mod.F90 \
          ngwf_cg_mod.F90 ngwf_data_mod.F90 ngwf_diis_mod.F90 \
          ngwf_gradient_mod.F90 ngwf_rep_mod.F90 ngwfs_mod.F90 \
          nho_mod.F90 npa_mod.F90 norm_conserv_pseudo_mod.F90 \
          nonsc_forces_mod.F90 onetep.F90 palser_mano_mod.F90 \
          paw_mod.F90 paw_shape_mod.F90 paw_xc_mod.F90 \
          parallel_strategy_mod.F90 pbc_corrections_mod.F90 \
          penalty_mod.F90 phonon_mod.F90 potential_mod.F90 \
          ppd_strategy_mod.F90 projectors_mod.F90 properties_mod.F90 \
          restart_mod.F90 rundat_blocks_mod.F90 rundat_mod.F90 \
          services_mod.F90 sparse_mod.F90 sparse_init_mod.F90 \
          spherical_wave_mod.F90 tddft_mod.F90 timer_mod.F90 \
          transport_mod.F90 tssearch_algor_mod.F90 tssearch_mod.F90 \
          tssearch_utils_mod.F90 utils_mod.F90 \
          vdwcorrection_mod.F90 visual_mod.F90 wrappers_mod.F90 \
          xc_component_mod.F90 xc_mod.F90 xc_funcs_mod.F90


OBJECTS = $(SOURCES:.F90=.o)

#
# General rules
#

.F90.o:
	$(F90) $(FFLAGS) $(EXTFLAGS) $(INCPATH) -c $<

#
# Specific rules
#

default:
	@echo 'No target specified!'

onetep: $(OBJECTS)
	$(F90) $(FFLAGS) $(EXTFLAGS) -o $(EXECDIR)/onetep.$(ARCH) $(OBJECTS) $(LIBPATH) $(LIBS)

debug:: EXTFLAGS = -DDEBUG $(DEBUGFLAGS)
debug:: onetep

timing:: EXTFLAGS = $(EXTFLAGS) -DTIMINGS
timing:: onetep

profile:: EXTFLAGS = $(OPTFLAGS) $(PROFILEFLAGS)
profile:: onetep

warnings:: EXTFLAGS = $(OPTFLAGS) $(WARNINGFLAGS)
warnings:: onetep

#
# Dependencies
#

atom_mod.o: atom_mod.F90 comms_mod.o constants_mod.o ion_mod.o \
            norm_conserv_pseudo_mod.o paw_mod.o rundat_mod.o services_mod.o \
            spherical_wave_mod.o utils_mod.o xc_mod.o wrappers_mod.o

augmentation_mod.o: augmentation_mod.F90 cell_mod.o \
                    cell_grid_mod.o comms_mod.o constants_mod.o \
                    function_basis_mod.o fourier_mod.o geometry_mod.o \
                    ion_mod.o norm_conserv_pseudo_mod.o \
                    parallel_strategy_mod.o paw_mod.o projectors_mod.o \
                    rundat_mod.o sparse_mod.o timer_mod.o utils_mod.o xc_mod.o

bandstructure_mod.o: bandstructure_mod.F90 cell_mod.o comms_mod.o \
                     constants_mod.o function_basis_mod.o geometry_mod.o \
                     integrals_mod.o ion_mod.o ngwf_rep_mod.o \
                     parallel_strategy_mod.o projectors_mod.o \
                     norm_conserv_pseudo_mod.o rundat_mod.o sparse_mod.o \
                     utils_mod.o hubbard_build_mod.o hubbard_init_mod.o 

basis_mod.o: basis_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
             geometry_mod.o timer_mod.o utils_mod.o

chebyshev_rep_mod.o: chebyshev_rep_mod.F90 comms_mod.o constants_mod.o \
                     utils_mod.o

cdft_cg_mod.o: cdft_cg_mod.F90 cell_mod.o cell_grid_mod.o \
               comms_mod.o constants_mod.o \
               function_basis_mod.o hamiltonian_mod.o \
               hubbard_build_mod.o hubbard_init_mod.o ion_mod.o \
               kernel_mod.o ngwf_cg_mod.o ngwf_diis_mod.o ngwf_rep_mod.o  \
               palser_mano_mod.o parallel_strategy_mod.o restart_mod.o \
               rundat_mod.o services_mod.o sparse_mod.o timer_mod.o \
               utils_mod.o wrappers_mod.o

cell_mod.o: cell_mod.F90 comms_mod.o constants_mod.o geometry_mod.o ion_mod.o \
            rundat_mod.o utils_mod.o

cell_grid_mod.o: cell_grid_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
                 geometry_mod.o rundat_mod.o timer_mod.o utils_mod.o

classical_pot_mod.o:  classical_pot_mod.F90 cell_mod.o cell_grid_mod.o \
                      comms_mod.o constants_mod.o geometry_mod.o ion_mod.o \
                      rundat_mod.o services_mod.o timer_mod.o utils_mod.o

comms_mod.o: comms_mod.F90 constants_mod.o

conduction_mod.o: conduction_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o \
                  conduction_properties_mod.o constants_mod.o electronic_mod.o \
                  function_basis_mod.o hamiltonian_mod.o hubbard_build_mod.o \
                  ion_mod.o kernel_mod.o ngwf_cg_mod.o ngwfs_mod.o \
                  ngwf_rep_mod.o rundat_mod.o services_mod.o \
                  sparse_mod.o utils_mod.o

conduction_properties_mod.o: conduction_properties_mod.F90 cell_mod.o cell_grid_mod.o \
                  comms_mod.o constants_mod.o function_basis_mod.o \
                  hamiltonian_mod.o hubbard_build_mod.o ion_mod.o kernel_mod.o \
                  ngwfs_mod.o ngwf_rep_mod.o properties_mod.o \
                  rundat_mod.o services_mod.o sparse_mod.o utils_mod.o

constants_mod.o: constants_mod.F90

cutoff_coulomb_mod.o: cutoff_coulomb_mod.F90 cell_mod.o \
                      cell_grid_mod.o classical_pot_mod.o comms_mod.o \
                      constants_mod.o density_mod.o fourier_mod.o \
                      geometry_mod.o hartree_mod.o ion_mod.o \
                      norm_conserv_pseudo_mod.o parallel_strategy_mod.o \
                      paw_mod.o rundat_mod.o services_mod.o timer_mod.o \
                      utils_mod.o visual_mod.o


dense_mod.o: dense_mod.F90 comms_mod.o constants_mod.o rundat_mod.o \
             sparse_mod.o timer_mod.o utils_mod.o

density_mod.o: density_mod.F90 basis_mod.o cell_mod.o cell_grid_mod.o \
               comms_mod.o constants_mod.o fourier_mod.o function_basis_mod.o \
               geometry_mod.o integrals_mod.o ion_mod.o parallel_strategy_mod.o \
               rundat_mod.o services_mod.o sparse_mod.o \
               timer_mod.o utils_mod.o visual_mod.o

electronic_mod.o: electronic_mod.F90 augmentation_mod.o cell_grid_mod.o \
                  cell_mod.o classical_pot_mod.o comms_mod.o constants_mod.o \
                  cutoff_coulomb_mod.o density_mod.o fourier_mod.o \
                  function_basis_mod.o hamiltonian_mod.o hartree_mod.o \
                  hubbard_build_mod.o integrals_mod.o ion_mod.o kernel_mod.o \
                  kernel_diis_mod.o lnv_mod.o ngwf_rep_mod.o \
                  norm_conserv_pseudo_mod.o palser_mano_mod.o paw_mod.o \
                  penalty_mod.o potential_mod.o restart_mod.o rundat_mod.o \
                  services_mod.o sparse_mod.o utils_mod.o xc_mod.o

eigenstates_mod.o: eigenstates_mod.F90 augmentation_mod.o basis_mod.o \
                   cell_mod.o cell_grid_mod.o comms_mod.o \
                   constants_mod.o dense_mod.o density_mod.o \
                   function_basis_mod.o integrals_mod.o ion_mod.o \
                   ngwf_rep_mod.o optics_mod.o parallel_strategy_mod.o \
                   rundat_mod.o sparse_mod.o timer_mod.o utils_mod.o \
                   visual_mod.o

energy_and_force_mod.o: energy_and_force_mod.F90 augmentation_mod.o \
                        cdft_cg_mod.o cell_grid_mod.o cell_mod.o \
                        classical_pot_mod.o comms_mod.o conduction_mod.o \
                        constants_mod.o cutoff_coulomb_mod.o \
                        dense_mod.o density_mod.o electronic_mod.o ewald_mod.o \
                        forces_mod.o fourier_mod.o function_basis_mod.o \
                        geometry_mod.o hf_exchange_mod.o hubbard_build_mod.o \
                        hubbard_init_mod.o \
                        ion_mod.o is_solvation_mod.o kernel_mod.o \
                        multigrid_methods_mod.o ngwf_cg_mod.o ngwf_diis_mod.o \
                        ngwf_rep_mod.o ngwfs_mod.o norm_conserv_pseudo_mod.o \
                        parallel_strategy_mod.o paw_mod.o \
                        pbc_corrections_mod.o potential_mod.o \
                        ppd_strategy_mod.o projectors_mod.o properties_mod.o \
                        rundat_mod.o services_mod.o sparse_mod.o tddft_mod.o \
                        utils_mod.o vdwcorrection_mod.o xc_mod.o

enthalpy_mod.o: enthalpy_mod.F90 cell_grid_mod.o cell_mod.o comms_mod.o \
                constants_mod.o integrals_mod.o rundat_mod.o utils_mod.o

esdf_key_mod.o: esdf_key_mod.F90

esdf_mod.o: esdf_mod.F90 constants_mod.o esdf_key_mod.o comms_mod.o utils_mod.o

ewald_mod.o: ewald_mod.F90 cell_mod.o classical_pot_mod.o comms_mod.o constants_mod.o \
             ion_mod.o rundat_mod.o timer_mod.o utils_mod.o

finite_differences_mod.o: finite_differences_mod.F90 cell_grid_mod.o \
             comms_mod.o constants_mod.o geometry_mod.o timer_mod.o utils_mod.o

forces_mod.o: forces_mod.F90 augmentation_mod.o cell_mod.o cell_grid_mod.o \
              classical_pot_mod.o comms_mod.o constants_mod.o cutoff_coulomb_mod.o \
              density_mod.o ewald_mod.o function_basis_mod.o hamiltonian_mod.o \
              hubbard_build_mod.o ion_mod.o ngwf_rep_mod.o \
              norm_conserv_pseudo_mod.o paw_mod.o rundat_mod.o \
              sparse_mod.o timer_mod.o utils_mod.o vdwcorrection_mod.o

fourier_mod.o: cell_mod.o cell_grid_mod.o comms_mod.o constants_mod.o \
               rundat_mod.o timer_mod.o utils_mod.o

function_basis_mod.o: basis_mod.o cell_mod.o comms_mod.o constants_mod.o \
                      fourier_mod.o geometry_mod.o hubbard_init_mod.o \
                      parallel_strategy_mod.o ion_mod.o rundat_mod.o \
                      sparse_mod.o spherical_wave_mod.o timer_mod.o \
                      utils_mod.o

gaunt_coeff_mod.o: gaunt_coeff_mod.F90 cell_mod.o constants_mod.o utils_mod.o

geometry_mod.o: geometry_mod.F90 constants_mod.o

geometry_deloc_algor.o: geometry_deloc_algor.F90 geometry_deloc_utils.o

geometry_deloc_utils.o: geometry_deloc_utils.F90 cell_mod.o comms_mod.o \
                        constants_mod.o ion_mod.o rundat_mod.o services_mod.o

geometry_optimiser_mod.o: geometry_optimiser_mod.F90 cell_mod.o classical_pot_mod.o \
                          comms_mod.o constants_mod.o energy_and_force_mod.o forces_mod.o \
                          geometry_mod.o geometry_deloc_algor.o \
                          geometry_deloc_utils.o ion_mod.o rundat_mod.o \
                          services_mod.o utils_mod.o wrappers_mod.o

hamiltonian_mod.o: hamiltonian_mod.F90 eigenstates_mod.o augmentation_mod.o cell_mod.o \
                   cell_grid_mod.o comms_mod.o constants_mod.o \
                   cutoff_coulomb_mod.o dense_mod.o density_mod.o \
                   enthalpy_mod.o function_basis_mod.o hartree_mod.o \
                   hf_exchange_mod.o hubbard_build_mod.o integrals_mod.o \
                   ion_mod.o is_smeared_ions_mod.o ngwf_rep_mod.o \
                   norm_conserv_pseudo_mod.o paw_mod.o projectors_mod.o \
                   rundat_mod.o sparse_mod.o timer_mod.o utils_mod.o \
                   vdwcorrection_mod.o wrappers_mod.o xc_component_mod.o xc_mod.o

hartree_mod.o: hartree_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o \
               constants_mod.o fourier_mod.o geometry_mod.o \
               is_smeared_ions_mod.o is_solvation_mod.o \
               multigrid_methods_mod.o \
               pbc_corrections_mod.o rundat_mod.o timer_mod.o utils_mod.o

hf_exchange_mod.o: hf_exchange_mod.F90 basis_mod.o cell_mod.o \
                   chebyshev_rep_mod.o constants_mod.o \
                   comms_mod.o fourier_mod.o function_basis_mod.o \
                   geometry_mod.o ion_mod.o parallel_strategy_mod.o \
                   rundat_mod.o sparse_mod.o spherical_wave_mod.o timer_mod.o \
                   utils_mod.o wrappers_mod.o xc_mod.o

hubbard_build_mod.o: hubbard_build_mod.F90 optics_mod.o \
                     augmentation_mod.o \
                     cell_mod.o comms_mod.o \
                     constants_mod.o fourier_mod.o function_basis_mod.o \
                     geometry_mod.o hubbard_init_mod.o integrals_mod.o \
                     ion_mod.o parallel_strategy_mod.o paw_mod.o \
                     projectors_mod.o restart_mod.o rundat_mod.o \
                     sparse_mod.o timer_mod.o utils_mod.o wrappers_mod.o density_mod.o  

hubbard_init_mod.o: hubbard_init_mod.F90 cell_mod.o comms_mod.o \
                    constants_mod.o rundat_mod.o utils_mod.o

integrals_mod.o: integrals_mod.F90 basis_mod.o cell_mod.o cell_grid_mod.o \
                 comms_mod.o constants_mod.o fourier_mod.o geometry_mod.o \
                 kinetic_mod.o potential_mod.o rundat_mod.o sparse_mod.o \
                 timer_mod.o utils_mod.o function_basis_mod.o

ion_mod.o: ion_mod.F90 constants_mod.o geometry_mod.o


is_poisson_mod.o: is_poisson_mod.F90 cell_grid_mod.o \
                  constants_mod.o \
                  utils_mod.o rundat_mod.o visual_mod.o

is_smeared_ions_mod.o: is_smeared_ions_mod.F90 cell_mod.o cell_grid_mod.o \
                       comms_mod.o constants_mod.o geometry_mod.o ion_mod.o \
                       integrals_mod.o multigrid_methods_mod.o \
                       parallel_strategy_mod.o rundat_mod.o utils_mod.o

is_solvation_mod.o: is_solvation_mod.F90 cell_grid_mod.o cell_mod.o \
                    comms_mod.o constants_mod.o finite_differences_mod.o \
                    fourier_mod.o geometry_mod.o \
                    integrals_mod.o multigrid_methods_mod.o is_poisson_mod.o \
                    is_smeared_ions_mod.o rundat_mod.o utils_mod.o

kernel_mod.o: kernel_mod.F90 augmentation_mod.o cell_mod.o comms_mod.o \
              constants_mod.o dense_mod.o function_basis_mod.o integrals_mod.o \
              rundat_mod.o sparse_mod.o timer_mod.o utils_mod.o

kernel_diis_mod.o: kernel_diis_mod.F90 cell_grid_mod.o comms_mod.o \
                   constants_mod.o dense_mod.o function_basis_mod.o \
                   hamiltonian_mod.o hubbard_build_mod.o ion_mod.o \
                   kernel_mod.o ngwf_rep_mod.o rundat_mod.o cell_mod.o \
                   sparse_mod.o timer_mod.o utils_mod.o vdwcorrection_mod.o

kinetic_mod.o: kinetic_mod.F90 basis_mod.o cell_mod.o \
               constants_mod.o fourier_mod.o function_basis_mod.o \
               utils_mod.o

lnv_mod.o: lnv_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o constants_mod.o \
           function_basis_mod.o hamiltonian_mod.o hubbard_build_mod.o \
           ion_mod.o kernel_mod.o ngwf_rep_mod.o restart_mod.o \
           rundat_mod.o services_mod.o sparse_mod.o timer_mod.o utils_mod.o

md_mod.o: md_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
          energy_and_force_mod.o ion_mod.o md_ions_mod.o \
          md_thermostat_mod.o restart_mod.o \
          rundat_mod.o utils_mod.o services_mod.o

md_ions_mod.o: md_ions_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
          ion_mod.o rundat_mod.o utils_mod.o services_mod.o

md_thermostat_mod.o: md_thermostat_mod.F90 cell_mod.o comms_mod.o \
          constants_mod.o md_ions_mod.o rundat_mod.o utils_mod.o \
          services_mod.o


multigrid_mod.o: comms_mod.o constants_mod.o timer_mod.o utils_mod.o

multigrid_methods_mod.o: multigrid_methods_mod.F90 cell_mod.o cell_grid_mod.o \
                     comms_mod.o constants_mod.o finite_differences_mod.o \
                     geometry_mod.o integrals_mod.o multigrid_mod.o \
                     rundat_mod.o timer_mod.o utils_mod.o \
                     visual_mod.o

ngwf_cg_mod.o: ngwf_cg_mod.F90 eigenstates_mod.o augmentation_mod.o cell_mod.o cell_grid_mod.o \
               comms_mod.o constants_mod.o electronic_mod.o forces_mod.o \
               function_basis_mod.o geometry_mod.o hamiltonian_mod.o \
               hubbard_build_mod.o ion_mod.o integrals_mod.o \
               kernel_mod.o ngwf_gradient_mod.o \
               ngwf_rep_mod.o nho_mod.o nonsc_forces_mod.o \
               norm_conserv_pseudo_mod.o projectors_mod.o restart_mod.o \
               rundat_mod.o services_mod.o sparse_mod.o timer_mod.o \
               utils_mod.o visual_mod.o wrappers_mod.o

ngwf_diis_mod.o: ngwf_diis_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o \
               constants_mod.o electronic_mod.o function_basis_mod.o \
               geometry_mod.o hamiltonian_mod.o hubbard_build_mod.o ion_mod.o \
               integrals_mod.o ngwf_gradient_mod.o ngwf_rep_mod.o nho_mod.o \
               norm_conserv_pseudo_mod.o \
               projectors_mod.o restart_mod.o \
               rundat_mod.o services_mod.o sparse_mod.o timer_mod.o \
               utils_mod.o visual_mod.o wrappers_mod.o

ngwf_data_mod.o: ngwf_data_mod.F90 constants_mod.o

ngwf_gradient_mod.o: ngwf_gradient_mod.F90 augmentation_mod.o basis_mod.o \
                     cell_mod.o cell_grid_mod.o comms_mod.o constants_mod.o \
                     fourier_mod.o function_basis_mod.o \
                     geometry_mod.o hf_exchange_mod.o hubbard_build_mod.o \
                     ion_mod.o kernel_mod.o kernel_diis_mod.o \
                     kinetic_mod.o ngwf_rep_mod.o norm_conserv_pseudo_mod.o \
                     paw_mod.o potential_mod.o projectors_mod.o \
                     rundat_mod.o sparse_mod.o timer_mod.o utils_mod.o

ngwf_rep_mod.o: cell_mod.o comms_mod.o constants_mod.o ion_mod.o rundat_mod.o \
                sparse_init_mod.o sparse_mod.o utils_mod.o

ngwfs_mod.o: ngwfs_mod.F90 atom_mod.o basis_mod.o cell_mod.o comms_mod.o \
             constants_mod.o density_mod.o function_basis_mod.o \
             geometry_mod.o ion_mod.o ngwf_data_mod.o parallel_strategy_mod.o \
             projectors_mod.o restart_mod.o rundat_mod.o services_mod.o \
             spherical_wave_mod.o utils_mod.o

nho_mod.o: nho_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
           function_basis_mod.o geometry_mod.o ion_mod.o \
           parallel_strategy_mod.o rundat_mod.o sparse_mod.o timer_mod.o \
           utils_mod.o wrappers_mod.o

nonsc_forces_mod.o: nonsc_forces_mod.F90 cell_mod.o comms_mod.o \
             constants_mod.o function_basis_mod.o \
             integrals_mod.o kinetic_mod.o \
             parallel_strategy_mod.o sparse_mod.o timer_mod.o utils_mod.o

norm_conserv_pseudo_mod.o: norm_conserv_pseudo_mod.F90 cell_mod.o \
                           cell_grid_mod.o classical_pot_mod.o comms_mod.o \
                           constants_mod.o fourier_mod.o function_basis_mod.o \
                           geometry_mod.o ion_mod.o is_smeared_ions_mod.o \
                           parallel_strategy_mod.o pbc_corrections_mod.o \
                           projectors_mod.o rundat_mod.o services_mod.o \
                           sparse_mod.o spherical_wave_mod.o timer_mod.o \
                           utils_mod.o xc_mod.o

npa_mod.o: npa_mod.F90 cell_mod.o comms_mod.o constants_mod.o dense_mod.o \
           function_basis_mod.o ion_mod.o kernel_mod.o ngwf_data_mod.o \
           ngwf_rep_mod.o parallel_strategy_mod.o  rundat_mod.o \
           sparse_mod.o timer_mod.o utils_mod.o wrappers_mod.o

optics_mod.o: optics_mod.F90 augmentation_mod.o basis_mod.o comms_mod.o \
              constants_mod.o dense_mod.o function_basis_mod.o geometry_mod.o \
              integrals_mod.o ion_mod.o ngwf_rep_mod.o \
              norm_conserv_pseudo_mod.o rundat_mod.o cell_mod.o sparse_mod.o \
              timer_mod.o utils_mod.o

onetep.o: onetep.F90 comms_mod.o constants_mod.o energy_and_force_mod.o \
          esdf_mod.o ewald_mod.o geometry_optimiser_mod.o hubbard_init_mod.o \
          ion_mod.o md_mod.o phonon_mod.o \
          rundat_mod.o rundat_blocks_mod.o \
          services_mod.o spherical_wave_mod.o timer_mod.o tssearch_mod.o \
          utils_mod.o vdwcorrection_mod.o

palser_mano_mod.o: palser_mano_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
                   rundat_mod.o services_mod.o sparse_mod.o timer_mod.o

parallel_strategy_mod.o: parallel_strategy_mod.F90 cell_mod.o comms_mod.o \
                         constants_mod.o geometry_mod.o hubbard_init_mod.o \
                         ion_mod.o rundat_mod.o utils_mod.o

paw_mod.o: paw_mod.F90 basis_mod.o classical_pot_mod.o cell_mod.o \
           cell_grid_mod.o constants_mod.o comms_mod.o ion_mod.o fourier_mod.o \
           function_basis_mod.o gaunt_coeff_mod.o geometry_mod.o \
           paw_shape_mod.o paw_xc_mod.o parallel_strategy_mod.o \
           projectors_mod.o rundat_mod.o services_mod.o sparse_mod.o \
           spherical_wave_mod.o timer_mod.o utils_mod.o xc_mod.o

paw_shape_mod.o: paw_shape_mod.F90 basis_mod.o cell_mod.o cell_grid_mod.o \
                 constants_mod.o fourier_mod.o geometry_mod.o services_mod.o \
                 spherical_wave_mod.o utils_mod.o

paw_xc_mod.o: paw_xc_mod.F90 constants_mod.o gaunt_coeff_mod.o \
              services_mod.o xc_mod.o

pbc_corrections_mod.o: pbc_corrections_mod.F90 cell_mod.o cell_grid_mod.o \
                       comms_mod.o constants_mod.o fourier_mod.o \
                       geometry_mod.o ion_mod.o rundat_mod.o timer_mod.o \
                       utils_mod.o

penalty_mod.o: penalty_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o \
               constants_mod.o function_basis_mod.o hamiltonian_mod.o \
               hubbard_build_mod.o integrals_mod.o ion_mod.o kernel_mod.o \
               ngwf_rep_mod.o restart_mod.o rundat_mod.o services_mod.o \
               sparse_mod.o timer_mod.o utils_mod.o vdwcorrection_mod.o

phonon_mod.o: phonon_mod.F90 dense_mod.o constants_mod.o energy_and_force_mod.o \
              ion_mod.o comms_mod.o cell_mod.o utils_mod.o rundat_mod.o

potential_mod.o: potential_mod.F90 basis_mod.o cell_mod.o cell_grid_mod.o \
                 comms_mod.o constants_mod.o fourier_mod.o \
                 function_basis_mod.o parallel_strategy_mod.o rundat_mod.o \
                 timer_mod.o utils_mod.o

ppd_strategy_mod.o: ppd_strategy_mod.F90 cell_mod.o comms_mod.o \
                    constants_mod.o fourier_mod.o geometry_mod.o \
                    ion_mod.o rundat_mod.o

projectors_mod.o: projectors_mod.F90 basis_mod.o cell_mod.o comms_mod.o \
                  constants_mod.o fourier_mod.o function_basis_mod.o \
                  geometry_mod.o integrals_mod.o \
                  ion_mod.o parallel_strategy_mod.o \
                  rundat_mod.o services_mod.o sparse_mod.o \
                  spherical_wave_mod.o timer_mod.o utils_mod.o

properties_mod.o: properties_mod.F90 augmentation_mod.o bandstructure_mod.o \
                  basis_mod.o cell_mod.o cell_grid_mod.o comms_mod.o \
                  constants_mod.o dense_mod.o density_mod.o eigenstates_mod.o \
                  function_basis_mod.o geometry_mod.o hamiltonian_mod.o \
                  hubbard_build_mod.o integrals_mod.o ion_mod.o kernel_mod.o \
                  md_ions_mod.o ngwf_rep_mod.o  npa_mod.o \
                  parallel_strategy_mod.o paw_mod.o rundat_mod.o sparse_mod.o \
                  spherical_wave_mod.o timer_mod.o transport_mod.o \
                  utils_mod.o visual_mod.o wrappers_mod.o

restart_mod.o: restart_mod.F90 basis_mod.o cell_mod.o comms_mod.o \
               constants_mod.o dense_mod.o fourier_mod.o \
               function_basis_mod.o geometry_mod.o integrals_mod.o \
               ion_mod.o kernel_mod.o ngwf_rep_mod.o \
               norm_conserv_pseudo_mod.o parallel_strategy_mod.o \
               projectors_mod.o rundat_mod.o services_mod.o \
               spherical_wave_mod.o sparse_mod.o timer_mod.o \
               utils_mod.o 

rundat_mod.o: rundat_mod.F90 comms_mod.o constants_mod.o esdf_mod.o utils_mod.o

rundat_blocks_mod.o: rundat_blocks_mod.F90 classical_pot_mod.o \
                     comms_mod.o constants_mod.o esdf_mod.o geometry_mod.o \
                     hubbard_init_mod.o ion_mod.o md_thermostat_mod.o \
                     ppd_strategy_mod.o rundat_mod.o services_mod.o \
                     cell_mod.o utils_mod.o transport_mod.o vdwcorrection_mod.o

services_mod.o: services_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
                ion_mod.o rundat_mod.o utils_mod.o

sparse_mod.o: sparse_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
              ion_mod.o parallel_strategy_mod.o rundat_mod.o timer_mod.o \
              utils_mod.o

sparse_init_mod.o: sparse_init_mod.F90 cell_mod.o constants_mod.o comms_mod.o \
                   ion_mod.o parallel_strategy_mod.o rundat_mod.o \
                   sparse_mod.o timer_mod.o utils_mod.o

spherical_wave_mod.o: spherical_wave_mod.F90 cell_mod.o \
                      comms_mod.o constants_mod.o fourier_mod.o \
                      geometry_mod.o timer_mod.o utils_mod.o

tddft_mod.o: tddft_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o \
             constants_mod.o function_basis_mod.o \
             geometry_mod.o hamiltonian_mod.o integrals_mod.o ion_mod.o \
             ngwf_rep_mod.o properties_mod.o \
             rundat_mod.o sparse_mod.o timer_mod.o \
             utils_mod.o vdwcorrection_mod.o wrappers_mod.o xc_mod.o

timer_mod.o: timer_mod.F90 comms_mod.o constants_mod.o rundat_mod.o

transport_mod.o: transport_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
                 dense_mod.o function_basis_mod.o ion_mod.o ngwf_rep_mod.o \
                 parallel_strategy_mod.o rundat_mod.o sparse_mod.o \
                 timer_mod.o utils_mod.o  

tssearch_algor_mod.o: tssearch_algor_mod.F90 tssearch_utils_mod.o

tssearch_mod.o: tssearch_mod.F90 cell_mod.o comms_mod.o constants_mod.o \
                esdf_mod.o ion_mod.o rundat_mod.o rundat_blocks_mod.o \
                tssearch_algor_mod.o tssearch_utils_mod.o

tssearch_utils_mod.o: tssearch_utils_mod.F90 cell_mod.o comms_mod.o \
                      constants_mod.o energy_and_force_mod.o forces_mod.o \
                      ion_mod.o utils_mod.o

utils_mod.o: utils_mod.F90 comms_mod.o constants_mod.o

vdwcorrection_mod.o: vdwcorrection_mod.F90 cell_mod.o comms_mod.o \
                     constants_mod.o geometry_mod.o ion_mod.o rundat_mod.o \
                     utils_mod.o

visual_mod.o: visual_mod.F90 basis_mod.o cell_mod.o cell_grid_mod.o \
              comms_mod.o constants_mod.o function_basis_mod.o geometry_mod.o \
              ion_mod.o parallel_strategy_mod.o rundat_mod.o utils_mod.o

wrappers_mod.o: wrappers_mod.F90 comms_mod.o constants_mod.o timer_mod.o \
                utils_mod.o

xc_component_mod.o: xc_component_mod.F90 augmentation_mod.o cell_mod.o \
		  cell_grid_mod.o comms_mod.o constants_mod.o \
          fourier_mod.o rundat_mod.o timer_mod.o utils_mod.o xc_funcs_mod.o

xc_funcs_mod.o: xc_funcs_mod.F90 constants_mod.o utils_mod.o

xc_mod.o: xc_mod.F90 cell_mod.o cell_grid_mod.o comms_mod.o constants_mod.o \
          fourier_mod.o rundat_mod.o timer_mod.o utils_mod.o xc_funcs_mod.o

