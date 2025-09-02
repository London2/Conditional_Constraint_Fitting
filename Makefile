conditional_constraint_example: conditional_constraint_example.cxx
	`root-config --cxx --cflags` -Isrc -o conditional_constraint_example conditional_constraint_example.cxx `root-config --libs` -Lsrc
