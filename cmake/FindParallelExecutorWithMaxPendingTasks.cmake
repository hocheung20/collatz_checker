find_path(ParallelExecutorWithMaxPendingTasks_INCLUDES
  NAMES
  ParallelExecutorWithMaxPendingTasks.hpp
  PATHS
  ${INCLUDE_INSTALL_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParallelExecutorWithMaxPendingTasks DEFAULT_MSG
                                  ParallelExecutorWithMaxPendingTasks_INCLUDES)
mark_as_advanced(ParallelExecutorWithMaxPendingTasks_INCLUDES)