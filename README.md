
**PhysiCell Version:**      1.6.0

Edited core/PhysiCell_phenotype{.h,.cpp} to add a flag (misnomer for test_b):
```
       bool execute_cycle_phase_entry_function;
```
use in core/PhysiCell_core.cpp, in divide():
```
	phenotype.execute_cycle_phase_entry_function = true; 
	child->phenotype.execute_cycle_phase_entry_function = true; 
```

Compile, run:
```
make test_b
test_b csc_prob_0.2.xml
test_b csc_prob_0.8.xml
```
  
