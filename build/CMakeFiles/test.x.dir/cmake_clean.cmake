FILE(REMOVE_RECURSE
  "CMakeFiles/test.x.dir/fortran/test.f90.o"
  "test.x.pdb"
  "test.x"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang Fortran)
  INCLUDE(CMakeFiles/test.x.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
