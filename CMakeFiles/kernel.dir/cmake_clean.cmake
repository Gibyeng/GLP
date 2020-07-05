file(REMOVE_RECURSE
  "libkernel.pdb"
  "libkernel.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/kernel.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
