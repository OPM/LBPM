===========================
Implementing a new LB model
===========================

While LBPM includes a range of fully-functioning lattice Boltzmann models, the commonly used
Bhatnager-Gross-Krook (BGK) model has been deliberately excluded. While the physical limitations
of this model are well-known, implementing the BGK model is an excellent way to understand
how to implement new LB models within the more general framework of LBPM. In this excercise
you will

* learn "what goes where"



* don't modify core data structures (unless you have a really good reason)

* re-use existing components whenever possible



