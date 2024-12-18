list(APPEND fortran_src
  sfs_psichi.f90
)

set(exe_name sfs_psichi.x)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(${exe_name} PRIVATE bacio::bacio_4
                                          w3emc::w3emc_4
					  sp::sp_4
                                          g2::g2_d
					  ip::ip_4)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
