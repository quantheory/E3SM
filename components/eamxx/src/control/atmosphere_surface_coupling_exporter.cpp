#include "atmosphere_surface_coupling_exporter.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <iomanip>

#include <array>

namespace scream
{
// =========================================================================================
SurfaceCouplingExporter::SurfaceCouplingExporter (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{

}
// =========================================================================================
void SurfaceCouplingExporter::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();       // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  const auto m2 = m*m;
  const auto s2 = s*s;
  auto Wm2 = W/m2;
  Wm2.set_string("W/m2");

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Qunit = kg/kg;
  Qunit.set_string("kg/kg");

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  FieldLayout scalar2d_layout     { {COL   },      {m_num_cols                 } };
  FieldLayout vector3d_layout     { {COL,CMP,LEV}, {m_num_cols, 2, m_num_levs  } };
  FieldLayout scalar3d_layout_mid { {COL,LEV},     {m_num_cols,    m_num_levs  } };
  FieldLayout scalar3d_layout_int { {COL,ILEV},    {m_num_cols,    m_num_levs+1} };

  constexpr int ps = Spack::n;

  // These fields are required for computation/exports
  add_field<Required>("p_int",                scalar3d_layout_int,  Pa,     grid_name);
  add_field<Required>("pseudo_density",       scalar3d_layout_mid,  Pa,     grid_name, ps);
  add_field<Required>("phis",                 scalar2d_layout,      m2/s2,  grid_name);
  add_field<Required>("p_mid",                scalar3d_layout_mid,  Pa,     grid_name, ps);
  add_field<Required>("qv",                   scalar3d_layout_mid,  Qunit,  grid_name, "tracers", ps);
  add_field<Required>("T_mid",                scalar3d_layout_mid,  K,      grid_name, ps);
  // TODO: Switch horiz_winds to using U and V, note right now there is an issue with when the subfields are created, so can't switch yet.
  add_field<Required>("horiz_winds",          vector3d_layout,      m/s,    grid_name);
  add_field<Required>("sfc_flux_dir_nir",     scalar2d_layout,      Wm2,    grid_name);
  add_field<Required>("sfc_flux_dir_vis",     scalar2d_layout,      Wm2,    grid_name);
  add_field<Required>("sfc_flux_dif_nir",     scalar2d_layout,      Wm2,    grid_name);
  add_field<Required>("sfc_flux_dif_vis",     scalar2d_layout,      Wm2,    grid_name);
  add_field<Required>("sfc_flux_sw_net" ,     scalar2d_layout,      Wm2,    grid_name);
  add_field<Required>("sfc_flux_lw_dn"  ,     scalar2d_layout,      Wm2,    grid_name);
  add_field<Required>("precip_liq_surf_mass", scalar2d_layout,      kg/m2,  grid_name);
  add_field<Required>("precip_ice_surf_mass", scalar2d_layout,      kg/m2,  grid_name);

  create_helper_field("Sa_z",       scalar2d_layout, grid_name);
  create_helper_field("Sa_u",       scalar2d_layout, grid_name); 
  create_helper_field("Sa_v",       scalar2d_layout, grid_name); 
  create_helper_field("Sa_tbot",    scalar2d_layout, grid_name); 
  create_helper_field("Sa_ptem",    scalar2d_layout, grid_name);
  create_helper_field("Sa_pbot",    scalar2d_layout, grid_name); 
  create_helper_field("Sa_shum",    scalar2d_layout, grid_name); 
  create_helper_field("Sa_dens",    scalar2d_layout, grid_name);
  create_helper_field("Sa_pslv",    scalar2d_layout, grid_name);
  create_helper_field("Faxa_rainl", scalar2d_layout, grid_name);
  create_helper_field("Faxa_snowl", scalar2d_layout, grid_name);
  create_helper_field("Faxa_swndr", scalar2d_layout, grid_name);
  create_helper_field("Faxa_swvdr", scalar2d_layout, grid_name);
  create_helper_field("Faxa_swndf", scalar2d_layout, grid_name);
  create_helper_field("Faxa_swvdf", scalar2d_layout, grid_name);
  create_helper_field("Faxa_swnet", scalar2d_layout, grid_name);
  create_helper_field("Faxa_lwdn",  scalar2d_layout, grid_name);
}
// =========================================================================================
void SurfaceCouplingExporter::create_helper_field (const std::string& name,
                                                   const FieldLayout& layout,
                                                   const std::string& grid_name)
{
  using namespace ekat::units;
  FieldIdentifier id(name,layout,Units::nondimensional(),grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation();
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}
// =========================================================================================
size_t SurfaceCouplingExporter::requested_buffer_size_in_bytes() const
{
  // Number of Reals needed by local views in the interface
  return Buffer::num_2d_vector_mid*m_num_cols*ekat::npack<Spack>(m_num_levs)*sizeof(Spack) +
         Buffer::num_2d_vector_int*m_num_cols*ekat::npack<Spack>(m_num_levs+1)*sizeof(Spack);
}
// =========================================================================================
void SurfaceCouplingExporter::init_buffers(const ATMBufferManager &buffer_manager)
{
  const int nlev_packs       = ekat::npack<Spack>(m_num_levs);
  const int nlevi_packs      = ekat::npack<Spack>(m_num_levs+1);

  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 2d views packed views
  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  m_buffer.dz = decltype(m_buffer.dz)(s_mem, m_num_cols, nlev_packs);
  s_mem += m_buffer.dz.size();
  m_buffer.z_mid = decltype(m_buffer.z_mid)(s_mem, m_num_cols, nlev_packs);
  s_mem += m_buffer.z_mid.size();
  m_buffer.z_int = decltype(m_buffer.z_int)(s_mem, m_num_cols, nlevi_packs);
  s_mem += m_buffer.z_int.size();

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);

  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SurfaceCouplingExporter.");
}
// =========================================================================================
void SurfaceCouplingExporter::setup_surface_coupling_data(const SCDataManager &sc_data_manager)
{
  m_num_cpl_exports    = sc_data_manager.get_num_cpl_fields();
  m_num_scream_exports = sc_data_manager.get_num_scream_fields();

  EKAT_ASSERT_MSG(m_num_scream_exports <= m_num_cpl_exports,
                  "Error! More SCREAM exports than actual cpl exports.\n");
  EKAT_ASSERT_MSG(m_num_cols == sc_data_manager.get_field_size(), "Error! Surface Coupling exports need to have size ncols.");

  // The export data is of size ncols,num_cpl_exports. All other data is of size num_scream_exports
  m_cpl_exports_view_h = decltype(m_cpl_exports_view_h) (sc_data_manager.get_field_data_ptr(),
                                                         m_num_cols, m_num_cpl_exports);
  m_cpl_exports_view_d = Kokkos::create_mirror_view(DefaultDevice(), m_cpl_exports_view_h);

  m_export_field_names = new name_t[m_num_scream_exports];
  std::memcpy(m_export_field_names, sc_data_manager.get_field_name_ptr(), m_num_scream_exports*32*sizeof(char));

  m_cpl_indices_view =
    decltype(m_cpl_indices_view)          (sc_data_manager.get_field_cpl_indices_ptr(),
                                           m_num_scream_exports);

  m_vector_components_view =
      decltype(m_vector_components_view)  (sc_data_manager.get_field_vector_components_ptr(),
                                           m_num_scream_exports);
  m_constant_multiple_view =
      decltype(m_constant_multiple_view)  (sc_data_manager.get_field_constant_multiple_ptr(),
                                           m_num_scream_exports);
  m_do_export_during_init_view = 
    decltype(m_do_export_during_init_view)(sc_data_manager.get_field_transfer_during_init_ptr(),
                                           m_num_scream_exports);

  m_column_info_d = decltype(m_column_info_d) ("m_info", m_num_scream_exports);
  m_column_info_h = Kokkos::create_mirror_view(m_column_info_d);
}
// =========================================================================================
void SurfaceCouplingExporter::initialize_impl (const RunType /* run_type */)
{
  bool any_initial_exports = false;

  for (int i=0; i<m_num_scream_exports; ++i) {

    std::string fname = m_export_field_names[i];
    EKAT_REQUIRE_MSG(has_helper_field(fname),"Error! Attempting to export "+fname+
                   " which has not been added as a helper field.\n");
    auto& field = m_helper_fields.at(fname);

    // Check that is valid
    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Export field view has not been allocated yet.\n");

    // Set view data ptr. Since the field could be only "Required", we
    // must use the unsafe version of the get_internal_view_data().
    m_column_info_h(i).data = field.get_internal_view_data_unsafe<Real>();

    // Get column info from field utility function
    get_col_info_for_surface_values(field.get_header_ptr(),
                                    m_vector_components_view(i),
                                    m_column_info_h(i).col_offset,
                                    m_column_info_h(i).col_stride);

    // Set constant multiple
    m_column_info_h(i).constant_multiple = m_constant_multiple_view(i);

    // Set whether or not this field should be exported during initialization
    m_column_info_h(i).transfer_during_initialization = m_do_export_during_init_view(i);
    if (m_do_export_during_init_view(i)) any_initial_exports = true;

    // Set index for referencing cpl data.
    m_column_info_h(i).cpl_indx = m_cpl_indices_view(i);
  }

  // Copy data to device for use in do_export()
  Kokkos::deep_copy(m_column_info_d, m_column_info_h);

  // Set the number of exports from eamxx or set to a constant, default type = EAMXX
  using vos_type = std::vector<std::string>;
  using vor_type = std::vector<Real>;
  m_export_source     = view_1d<DefaultDevice,ExportType>("",m_num_scream_exports);
  auto export_source_h = Kokkos::create_mirror_view(m_export_source);
  Kokkos::deep_copy(m_export_source,EAMXX);  // The default is that all export variables will be derived from the EAMxx state.
  m_num_eamxx_exports = m_num_scream_exports;
 
  if (m_params.isSublist("prescribed_constants")) {
    auto export_constant_params = m_params.sublist("prescribed_constants");
    EKAT_REQUIRE_MSG(export_constant_params.isParameter("fields"),"Error! surface_coupling_exporter::init - prescribed_constants does not have 'fields' parameter.");
    EKAT_REQUIRE_MSG(export_constant_params.isParameter("values"),"Error! surface_coupling_exporter::init - prescribed_constants does not have 'values' parameter.");
    auto export_constant_fields = export_constant_params.get<vos_type>("fields");
    auto export_constant_values = export_constant_params.get<vor_type>("values");
    EKAT_REQUIRE_MSG(export_constant_fields.size()==export_constant_values.size(),"Error! surface_coupling_exporter::init - prescribed_constants 'fields' and 'values' are not the same size");
    if (export_constant_fields.size()>0) {
      // Determine which fields need constants
      m_export_constants = view_1d<HostDevice,Real>("",m_num_scream_exports);
      for (int i=0; i<m_num_scream_exports; ++i) {  // TODO: This loop would probably be simpler if we just checked which "i" corresponded to each name in the fields list.
        std::string fname = m_export_field_names[i];
        auto loc = std::find(export_constant_fields.begin(),export_constant_fields.end(),fname);
        if (loc != export_constant_fields.end()) {
          const auto pos = loc-export_constant_fields.begin();
          export_source_h(i) = CONSTANT;
          m_num_const_exports += 1;
          m_num_eamxx_exports -= 1;
          m_export_constants(i) = export_constant_values[pos];
        }
      }
    }
  }
  // Copy host view back to device view
  Kokkos::deep_copy(m_export_source,export_source_h);
  // Final sanity check
  EKAT_REQUIRE_MSG(m_num_scream_exports = m_num_const_exports+m_num_eamxx_exports,"Error! surface_coupling_exporter - Something went wrong set the type of export for all variables.");
  EKAT_REQUIRE_MSG(m_num_eamxx_exports>=0,"Error! surface_coupling_exporter - The number of exports derived from EAMxx < 0, something must have gone wrong in assigning the types of exports for all variables.");

  // Perform initial export (if any are marked for export during initialization)
  if (any_initial_exports) do_export(0, true);
}
// =========================================================================================
void SurfaceCouplingExporter::run_impl (const double dt)
{
  do_export(dt);
}
// =========================================================================================
void SurfaceCouplingExporter::do_export(const double dt, const bool called_during_initialization)
{
  if (m_num_const_exports>0) {
    do_export_constant(dt,called_during_initialization);
  }
  if (m_num_eamxx_exports>0) {
    do_export_from_eamxx(dt,called_during_initialization);
  }

  // Finish up exporting vars
  do_export_to_cpl(called_during_initialization);
}
// =========================================================================================
void SurfaceCouplingExporter::do_export_constant(const double dt, const bool called_during_initialization)
{
  // Cycle through those fields that will be set to a constant value:
  auto export_source_h = Kokkos::create_mirror_view(m_export_source);
  Kokkos::deep_copy(export_source_h,m_export_source);
  for (int i=0; i<m_num_scream_exports; ++i) {
    if (export_source_h(i)==CONSTANT) {
      std::string fname = m_export_field_names[i];
      const auto field_view = m_helper_fields.at(fname).get_view<Real*>();
      Kokkos::deep_copy(field_view,m_export_constants(i));
    }
  }
  
}
// =========================================================================================
// This do_export handles all export variables that are derived from the EAMxx state.
// Important! This setup assumes the numerical order of export_cpl_indices as listed in
// /src/mct_coupling/scream_cpl_indices.F90
//
// If this order is changed or a new variable is added it is important to update the corresponding
// index query in the below.
void SurfaceCouplingExporter::do_export_from_eamxx(const double dt, const bool called_during_initialization)
{
  using PC = physics::Constants<Real>;

  const auto& p_int                = get_field_in("p_int").get_view<const Real**>();
  const auto& pseudo_density       = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& qv                   = get_field_in("qv").get_view<const Spack**>();
  const auto& T_mid                = get_field_in("T_mid").get_view<const Spack**>();
  // TODO: This will need to change if we ever switch from horiz_winds to U and V
  const auto& horiz_winds          = get_field_in("horiz_winds").get_view<const Real***>();
  const auto& p_mid                = get_field_in("p_mid").get_view<const Spack**>();
  const auto& phis                 = get_field_in("phis").get_view<const Real*>();
  const auto& sfc_flux_dir_nir     = get_field_in("sfc_flux_dir_nir").get_view<const Real*>();
  const auto& sfc_flux_dir_vis     = get_field_in("sfc_flux_dir_vis").get_view<const Real*>();
  const auto& sfc_flux_dif_nir     = get_field_in("sfc_flux_dif_nir").get_view<const Real*>();
  const auto& sfc_flux_dif_vis     = get_field_in("sfc_flux_dif_vis").get_view<const Real*>();
  const auto& sfc_flux_sw_net      = get_field_in("sfc_flux_sw_net" ).get_view<const Real*>();
  const auto& sfc_flux_lw_dn       = get_field_in("sfc_flux_lw_dn"  ).get_view<const Real*>();

  const auto& precip_liq_surf_mass = get_field_in("precip_liq_surf_mass").get_view<const Real*>();
  const auto& precip_ice_surf_mass = get_field_in("precip_ice_surf_mass").get_view<const Real*>();

  const auto Sa_z       = m_helper_fields.at("Sa_z").get_view<Real*>();
  const auto Sa_u       = m_helper_fields.at("Sa_u").get_view<Real*>();
  const auto Sa_v       = m_helper_fields.at("Sa_v").get_view<Real*>();
  const auto Sa_tbot    = m_helper_fields.at("Sa_tbot").get_view<Real*>();
  const auto Sa_ptem    = m_helper_fields.at("Sa_ptem").get_view<Real*>();
  const auto Sa_pbot    = m_helper_fields.at("Sa_pbot").get_view<Real*>();
  const auto Sa_shum    = m_helper_fields.at("Sa_shum").get_view<Real*>();
  const auto Sa_dens    = m_helper_fields.at("Sa_dens").get_view<Real*>();
  const auto Sa_pslv    = m_helper_fields.at("Sa_pslv").get_view<Real*>();
  const auto Faxa_rainl = m_helper_fields.at("Faxa_rainl").get_view<Real*>();
  const auto Faxa_snowl = m_helper_fields.at("Faxa_snowl").get_view<Real*>();
  const auto Faxa_swndr = m_helper_fields.at("Faxa_swndr").get_view<Real*>();
  const auto Faxa_swvdr = m_helper_fields.at("Faxa_swvdr").get_view<Real*>();
  const auto Faxa_swndf = m_helper_fields.at("Faxa_swndf").get_view<Real*>();
  const auto Faxa_swvdf = m_helper_fields.at("Faxa_swvdf").get_view<Real*>();
  const auto Faxa_swnet = m_helper_fields.at("Faxa_swnet").get_view<Real*>();
  const auto Faxa_lwdn  = m_helper_fields.at("Faxa_lwdn" ).get_view<Real*>();

  const auto dz    = m_buffer.dz;
  const auto z_int = m_buffer.z_int;
  const auto z_mid = m_buffer.z_mid;


  // Local copies, to deal with CUDA's handling of *this.
  const int  num_levs           = m_num_levs;
  const int  num_cols           = m_num_cols;

  // Preprocess exports
  const auto setup_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(num_cols, num_levs);
  Kokkos::parallel_for(setup_policy, KOKKOS_LAMBDA(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) {
    const int i = team.league_rank();

    const auto qv_i             = ekat::subview(qv, i);
    const auto u_wind_i         = ekat::subview(horiz_winds, i, 0); // TODO, when U and V work switch to using here.
    const auto v_wind_i         = ekat::subview(horiz_winds, i, 1);
    const auto T_mid_i          = ekat::subview(T_mid, i);
    const auto p_mid_i          = ekat::subview(p_mid, i);
    const auto p_int_i          = ekat::subview(p_int, i);
    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto dz_i             = ekat::subview(dz, i);
    const auto z_int_i          = ekat::subview(z_int, i);
    const auto z_mid_i          = ekat::subview(z_mid, i);

    // Compute vertical layer thickness
    PF::calculate_dz(team, pseudo_density_i, p_mid_i, T_mid_i, qv_i, dz_i);
    team.team_barrier();

    // Compute vertical layer heights (relative to ground surface rather than from sea level).
    // Use z_int(nlevs) = z_surf = 0.0.
    const Real z_surf = 0.0;
    PF::calculate_z_int(team, num_levs, dz_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, num_levs, z_int_i, z_mid_i);
    team.team_barrier();

    const auto s_qv_i = ekat::scalarize(qv_i);
    const auto s_dz_i = ekat::scalarize(dz_i);
    const auto s_z_mid_i = ekat::scalarize(z_mid_i);
    const auto s_pseudo_density_i = ekat::scalarize(pseudo_density_i);
    const auto s_p_mid_i = ekat::scalarize(p_mid_i);
    const auto s_T_mid_i = ekat::scalarize(T_mid_i);

    // Calculate air temperature at bottom of cell closest to the ground for PSL
    const Real T_int_bot = PF::calculate_surface_air_T(s_T_mid_i(num_levs-1),s_z_mid_i(num_levs-1));

    // Set the values in the helper fields which correspond to the exported variables
    if (m_export_source(0)==EAMXX) { Sa_z(i)    = s_z_mid_i(num_levs-1); }
    if (m_export_source(1)==EAMXX) { Sa_u(i)    = u_wind_i(num_levs-1); }
    if (m_export_source(2)==EAMXX) { Sa_v(i)    = v_wind_i(num_levs-1); }
    if (m_export_source(3)==EAMXX) { Sa_tbot(i) = s_T_mid_i(num_levs-1); }
    if (m_export_source(4)==EAMXX) { Sa_ptem(i) = PF::calculate_theta_from_T(s_T_mid_i(num_levs-1), s_p_mid_i(num_levs-1)); }
    if (m_export_source(5)==EAMXX) { Sa_pbot(i) = s_p_mid_i(num_levs-1);  }
    if (m_export_source(6)==EAMXX) { Sa_shum(i) = s_qv_i(num_levs-1); }
    if (m_export_source(7)==EAMXX) { Sa_dens(i) = PF::calculate_density(s_pseudo_density_i(num_levs-1), s_dz_i(num_levs-1)); }
    if (m_export_source(8)==EAMXX) { Sa_pslv(i) = PF::calculate_psl(T_int_bot, p_int_i(num_levs), phis(i)); }

    if (not called_during_initialization) {
      // Precipitation has units of kg/m2, and Faxa_rainl/snowl
      // need units mm/s. Here, 1000 converts m->mm, dt has units s, and
      // rho_h2o has units kg/m3.
      if (m_export_source(9)==EAMXX)  { Faxa_rainl(i) = precip_liq_surf_mass(i)/dt*(1000.0/PC::RHO_H2O); }
      if (m_export_source(10)==EAMXX) { Faxa_snowl(i) = precip_ice_surf_mass(i)/dt*(1000.0/PC::RHO_H2O); }
    }
    if (m_export_source(11)==EAMXX) { Faxa_swndr(i) = sfc_flux_dir_nir(i); }
    if (m_export_source(12)==EAMXX) { Faxa_swvdr(i) = sfc_flux_dir_vis(i); }
    if (m_export_source(13)==EAMXX) { Faxa_swndf(i) = sfc_flux_dif_nir(i); }
    if (m_export_source(14)==EAMXX) { Faxa_swvdf(i) = sfc_flux_dif_vis(i); }
    if (m_export_source(15)==EAMXX) { Faxa_swnet(i) = sfc_flux_sw_net(i); }
    if (m_export_source(16)==EAMXX) { Faxa_lwdn (i) = sfc_flux_lw_dn(i); }
  });
}
// =========================================================================================
void SurfaceCouplingExporter::do_export_to_cpl(const bool called_during_initialization)
{
  using policy_type = KT::RangePolicy;
  // Any field not exported by scream, or not exported
  // during initialization, is set to 0.0
  Kokkos::deep_copy(m_cpl_exports_view_d, 0.0);
  const auto cpl_exports_view_d = m_cpl_exports_view_d;
  const int  num_exports        = m_num_scream_exports;
  const int  num_cols           = m_num_cols;
  const auto col_info           = m_column_info_d;
  // Export to cpl data
  auto export_policy   = policy_type (0,num_exports*num_cols);
  Kokkos::parallel_for(export_policy, KOKKOS_LAMBDA(const int& i) {
    const int ifield = i / num_cols;
    const int icol   = i % num_cols;
    const auto& info = col_info(ifield);
    const auto offset = icol*info.col_stride + info.col_offset;

    // if this is during initialization, check whether or not the field should be exported
    bool do_export = (not called_during_initialization || info.transfer_during_initialization);
    if (do_export) {
      cpl_exports_view_d(icol,info.cpl_indx) = info.constant_multiple*info.data[offset];
    }
  });

  // Deep copy fields from device to cpl host array
  Kokkos::deep_copy(m_cpl_exports_view_h,m_cpl_exports_view_d);
}
// =========================================================================================
void SurfaceCouplingExporter::finalize_impl()
{

}
// =========================================================================================
} // namespace scream
