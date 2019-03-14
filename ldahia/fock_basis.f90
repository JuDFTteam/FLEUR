MODULE m_fock_basis
   USE m_juDFT

   CONTAINS
   !------------------------------------------------------------------------------
   !
   ! MODULE:  m_fock_basis
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module contains the subroutines to set up and perform operations
   !>  on the fock basis for the atomic hamiltonian used in DFT+Hubbard1    
   !>  The states are represented by bit strings where every bit determines 
   !>  wether a state is occupied or not                                    
   !>  The order of the states is as follows (labeled as [m,up/down]) where 
   !>  up/down is the spin index:                                          
   !>    ********************************************************************
   !>    [l,up], [l,down], [l-1,up], [l-1,down], .... , [-l,up], [-l,down]    
   !>    ******************************************************************** 
   !>  Here the last state [-l,down] sits on the bit with index 0           
   !>                                                                     
   !>  The goal is to keep much of the actual bit logic inside this         
   !>  module and work only with magnetic quantum numbers and spin indices  
   !>  outside of here                                           
   !
   ! REVISION HISTORY:
   ! February 2019 - Initial Version
   !------------------------------------------------------------------------------

   SUBROUTINE gen_fock_states(N,N_states,states)

      !finds all possible fock states with N electrons 
      !by generating the states in increasing order of the associated integer
      !The algorithm follows std::next_permutation in C++ but
      !adapted because we only work with bits

      IMPLICIT NONE

      INTEGER,             INTENT(IN)     :: N        !number of electrons
      INTEGER,             INTENT(IN)     :: N_states !number of one-particle states
      INTEGER,             INTENT(INOUT)  :: states(:)!Array of fock states

      LOGICAL next
      INTEGER length, i,j,k, N_f
      INTEGER count
      INTEGER state

      states(:) = 0

      IF(N.GT.N_states) CALL juDFT_error("Invalid occupation for DFT+Hubbard 1", &
                                        calledby="gen_fock_states",hint="This This is a BUG in FLEUR, please report")

      IF(N.EQ.0) THEN
         N_f = 0
         RETURN
      ENDIF

      state = 0
      DO i = 0,N-1
         state = ibset(state,i)
      ENDDO

      next = .true.

      count = 0
      DO WHILE(next)

         count = count+1
         states(count) = state

         !search for the first 01 in the bit string (read from the right side)
         j = 0
         DO WHILE(.NOT.(.NOT.btest(state,j+1).AND.btest(state,j)))
            j = j+1
            IF(j+1.EQ.N_states) next = .false.
         ENDDO

         !search for the first 1 to the right of j
         k = 0
         DO WHILE(k.LT.j.AND.(.NOT.btest(state,k)))
            k = k+1
         ENDDO

         !swap the 0 at j+1 with the 1 at k
         CALL swap(state,j+1,k)

         !reverse the order of the bit string from j to the end
         IF(j.GE.1) THEN
            DO i = 0, j/2
               CALL swap(state,i,j-i)
            END DO
         ENDIF

      ENDDO

      N_f = count

   END SUBROUTINE gen_fock_states

   SUBROUTINE find_state(state,array,N,N_states,ind)

      !finds the array index for a particular fock state
      !Uses the fact that the states are sorted in increasing order

      IMPLICIT NONE

      INTEGER,             INTENT(IN)  :: state
      INTEGER,             INTENT(IN)  :: array(:)
      INTEGER,             INTENT(IN)  :: N        !length of the array
      INTEGER,             INTENT(IN)  :: N_states !number of significant bits
      INTEGER,             INTENT(OUT) :: ind

      INTEGER i_bit, i
      INTEGER start_ind, end_ind
      LOGICAL bit

      start_ind = 1
      end_ind   = N 

      DO i_bit = N_states-1,0,-1 !loop over the significant bits

         bit = btest(state,i_bit)

         !search for the point where the significant bit switches
         i = start_ind
         DO WHILE (.NOT.btest(array(i),i_bit))
            i = i + 1
         ENDDO

         !restrict the search interval according 
         !to the bit of the state to be found
         IF(bit) THEN
            start_ind = i
         ELSE
            end_ind = i-1
         ENDIF

         !Check wether we found the state
         IF(start_ind.EQ.end_ind) THEN
            ind = start_ind
            EXIT
         END IF

      ENDDO

   END SUBROUTINE find_state

   SUBROUTINE find_transitions(ini_state,fin_state,N_states,N,allow_spin_flip,N_op,op,dist)

      !finds the possible transitions to get from ini_state to fin_state using
      !N pairs of creation and annihilation operators

      !allow_spin_flip determines wether the spin is allowed to change 
      !The output dist is given to avoid unnecessary double calling of the subroutines

      !TODO: clean the routine up (maybe generalize)

      IMPLICIT NONE

      INTEGER,                   INTENT(IN)  :: ini_state         !initial states
      INTEGER,                   INTENT(IN)  :: fin_state         !final state
      INTEGER,                   INTENT(IN)  :: N_states          !number of one-particle states
      INTEGER,                   INTENT(IN)  :: N                 !number of creation/annihilation operator pairs (only calculates for 1 and 2)
      LOGICAL,                   INTENT(IN)  :: allow_spin_flip   !are the operations allowed to change the spin

      INTEGER,                   INTENT(OUT) :: N_op              !number of possible operations
      INTEGER,ALLOCATABLE ,      INTENT(OUT) :: op(:,:)           !operations that allow this
      INTEGER,                   INTENT(OUT) :: dist              !difference between the two states

      INTEGER i, ann_ind, crea_ind, j, i_free
      INTEGER n_even, tmp,avail
      LOGICAL valid

      INTEGER, ALLOCATABLE :: diff(:),avail_bit(:)

      ALLOCATE(diff(N_states))
      ALLOCATE(avail_bit(N_states))

      !Look at the points where the two states differ:
      CALL distance(ini_state,fin_state,N_states,dist,diff)

      valid = .true.

      !We only look at atmost 2 pairs of creation/annihilation operatos
      IF(N.GE.3) THEN
         N_op = 0
         valid = .false.
      END IF

      !A transition is only possible if the states only differ in at most 2*N places
      !e.g. we can only change 4 indices with two pairs of creation and annihilation operators
      IF(dist.GT.2*N) THEN
         N_op = 0
         valid = .false.
      ELSE IF(.NOT.allow_spin_flip) THEN
         !If the operation is not allowed to change the spin look wether its still possible

         !spin indices are stored on odd/even bit indices
         !If changing the spin is not allowed the number of spin up/down particles cant change
         n_even = 0
         DO i = 1, dist
            IF(MOD(diff(i),2).EQ.0) THEN
               !Count the number of times a spin up is created(-1) annihilated(+1)
               IF(btest(ini_state,diff(i))) THEN
                  n_even = n_even+1
               ELSE
                  n_even = n_even-1
               END IF
            ENDIF
         ENDDO

         !This number must be 0
         IF(n_even.NE.0) THEN
            N_op = 0
            valid = .false.
         END IF

      END IF

      IF(valid) THEN

         !The maximum number of possible operations is maybe smaller. This should be large enouh in all cases
         !(just to be safe, it's not much space wasted) 

         ALLOCATE(op(2*N_states**N,2*N))

         op(:,:) = 0
         
         !Look for possible states if we have free pairs
         IF(dist.NE.2*N) THEN
            avail = 0
            avail_bit(:) = 0

            DO i = 0, N_states-1
               IF(btest(ini_state,i).AND.ALL(diff(:dist).NE.i)) THEN
                  avail = avail + 1
                  avail_bit(avail) = i
               END IF
            ENDDO
         ELSE
            avail = 1 !To make sure that we iterate once
         END IF
         
         N_op = 0

         DO i_free = 1,avail            

            N_op = N_op + 1

            ann_ind = N
            crea_ind = 0

            !Set the operations that are necessary to transform ini_state into fin_state
            IF(dist.NE.0) THEN
               DO i = 1, dist
                  IF(btest(ini_state,diff(i))) THEN
                     ann_ind = ann_ind + 1
                     op(N_op,ann_ind) = diff(i)
                  ELSE
                     crea_ind = crea_ind + 1
                     op(N_op,crea_ind) = diff(i)
                  END IF
               ENDDO
            END IF

            !Check if there are free pairs left
            IF(dist.NE.2*N) THEN
               !We have at most 2 free pairs in our case
               IF(N.EQ.2.AND.dist.EQ.0) THEN
                  DO j = i_free+1, avail
                     op(N_op,1) = avail_bit(j)
                     op(N_op,2) = avail_bit(i_free)
                     op(N_op,3) = avail_bit(i_free)
                     op(N_op,4) = avail_bit(j)
                     IF(j.NE.avail) N_op = N_op + 1
                  ENDDO
               ELSE 
                  !Otherwise we have only one free pair in every circumstance
                  op(N_op,crea_ind+1) = avail_bit(i_free)
                  op(N_op,ann_ind+1) = avail_bit(i_free)
               END IF
            END IF

         END DO

         !If we are dealing with two pairs we have permutations of the operation that are valid 
         IF(N.EQ.2) THEN
            IF(dist.EQ.0) N_op = N_op - 1 !In this case we had one addition to much
            tmp = N_op
            DO i = 1, tmp
               !Permute keeping the spin indices in the same pair
               N_op = N_op + 1 
               op(N_op,1) = op(i,2)
               op(N_op,2) = op(i,1)
               op(N_op,3) = op(i,4)
               op(N_op,4) = op(i,3)

               IF(MOD(op(i,3),2).EQ.MOD(op(i,4),2).OR.allow_spin_flip) THEN
                  !Permute mixing the spin indices (if they are equal) 
                  !The switch allow_spin_flip is here to keep it consistent, but it shouldnt be needed in our case
                  N_op = N_op + 1 
                  op(N_op,1) = op(i,1)
                  op(N_op,2) = op(i,2)
                  op(N_op,3) = op(i,4)
                  op(N_op,4) = op(i,3)

                  N_op = N_op + 1
                  op(N_op,1) = op(i,2)
                  op(N_op,2) = op(i,1)
                  op(N_op,3) = op(i,3)
                  op(N_op,4) = op(i,4)
               END IF
            ENDDO
         END IF
      END IF

   END SUBROUTINE

   SUBROUTINE bit_to_op(ind,n,l,m,spin)

      !Returns the corresponding m and spin index for a bit index

      INTEGER,    INTENT(IN)  :: ind(:)
      INTEGER,    INTENT(IN)  :: n
      INTEGER,    INTENT(IN)  :: l
      INTEGER,    INTENT(OUT) :: m(:)
      INTEGER,    INTENT(OUT) :: spin(:)

      INTEGER i
      DO i = 1, n 
         m(i) = INT(ind(i)/2)-l
         spin(i) = MOD(ind(i),2)
      ENDDO

   END SUBROUTINE


   !
   !OPERATIONS ON A FOCK-STATE OR A EIGENSTATE:
   !
   SUBROUTINE c_mu(dag,m,spin,l,N_states,N_in,N_out,vec_in,basis_in,vec_out,basis_out)

      !This module calculates the vector resulting from applying a creation/annihilation operator
      !to a vector in the fock-basis

      !Because the vectors have different lengths they are separated into two arguments
      !The basis_out has to be given and is the fock basis for N +- 1 electrons

      IMPLICIT NONE

      LOGICAL,             INTENT(IN)     :: dag
      INTEGER,             INTENT(IN)     :: m
      INTEGER,             INTENT(IN)     :: spin
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: N_in
      INTEGER,             INTENT(IN)     :: N_out
      INTEGER,             INTENT(IN)     :: N_states

      REAL,                INTENT(IN)     :: vec_in(:)
      REAL, ALLOCATABLE,   INTENT(OUT)    :: vec_out(:)

      INTEGER,             INTENT(IN)     :: basis_in(:)
      INTEGER,             INTENT(IN)     :: basis_out(:)

      INTEGER i, state, ind

      ALLOCATE(vec_out(N_out))
      vec_out(:) = 0.0
      !
      !Apply the creation/annihilation operator to every basis state and find it in basis_out
      !Then set the corresponding value of vec_in at this position in vec_out
      !
      DO i = 1, N_in
         CALL c(dag,m,spin,l,basis_in(i),state)
         !
         ! Only proceed if we do not get the vacuum state
         !
         IF(state.NE.0) THEN
            !
            !Find the new state in basis_out
            !
            CALL find_state(state,basis_out,N_out,N_states,ind)
            !
            !write the value of the old vector into the new vector
            !
            vec_out(ind) = vec_in(i)
         END IF
      ENDDO
   END SUBROUTINE c_mu

   SUBROUTINE c(dag,m,spin,l,state_in,state_out)

      !This is the creation/annihilation operator for given m and spin
      !Which of the two it is, is determined by dag (true for creation)

      !This is only used to calculate the <nu|c_[m,s]|mu> matrix elements
      !in the lehmann representation of the interacting green's function

      LOGICAL,      INTENT(IN)   :: dag
      INTEGER,      INTENT(IN)   :: m
      INTEGER,      INTENT(IN)   :: spin
      INTEGER,      INTENT(IN)   :: l
      INTEGER,      INTENT(IN)   :: state_in   
      INTEGER,      INTENT(OUT)  :: state_out 

      INTEGER i_bit

      !Find the corresponding bit index for [m,spin]
      i_bit = 2*(m+l)+ spin

      IF((btest(state_in,i_bit).AND.dag).OR.(.NOT.btest(state_in,i_bit).AND..NOT.dag)) THEN
         !Set the state to the vacuum state if the operation is not allowed
         state = 0
      ELSE
         IF(dag) THEN
            state_out = ibset(state_in,i_bit)
         ELSE
            state_out = ibclr(state_in,i_bit)
         END IF
      ENDIF

   END SUBROUTINE c
   !
   !ELEMENTAL BIT OPERATIONS ON THE STATE:
   !
   SUBROUTINE distance(N1,N2,states,diff,diff_bit)

      !determine the number of different indices between two states

      IMPLICIT NONE

      INTEGER,       INTENT(IN)  :: N1       !state 1
      INTEGER,       INTENT(IN)  :: N2       !state 2
      INTEGER,       INTENT(IN)  :: states   !number of bits/one-particle states
      INTEGER,       INTENT(OUT) :: diff     !number of different bits 
      INTEGER,       INTENT(OUT) :: diff_bit(states) !index of these bits
      INTEGER i

      diff = 0

      diff_bit(:) = 0 

      DO i = 0, states-1

         IF(XOR(btest(N1,i),btest(N2,i))) THEN
            diff = diff + 1
            diff_bit(diff) = i
         END IF

      ENDDO

   END SUBROUTINE distance

   SUBROUTINE swap(state,i1,i2)

      IMPLICIT NONE

      INTEGER,          INTENT(INOUT) ::  state
      INTEGER,          INTENT(IN) :: i1,i2

      LOGICAL tmp

      tmp = btest(state,i1)
      IF(btest(state,i2)) THEN
         state = ibset(state,i1)
      ELSE
         state = ibclr(state,i1)
      ENDIF
      IF(tmp) THEN
         state = ibset(state,i2)
      ELSE
         state = ibclr(state,i2)
      ENDIF

   END SUBROUTINE swap
   !
   !SOME IO-ROUTINES (JUST FOR DEBUG):
   !
   SUBROUTINE write_state(state,N_states,str)

      IMPLICIT NONE

      INTEGER,                 INTENT(IN)  :: state
      INTEGER,                 INTENT(IN)  :: N_states
      CHARACTER(len=N_states), INTENT(OUT):: str
      INTEGER i

      DO i = 1, N_states
         IF(btest(state,N_states-i)) THEN
            str(i:i) = "1"
         ELSE
            str(i:i) = "0"
         END IF
      ENDDO

   END SUBROUTINE write_state
   !
   !MATHEMATICAL FUNCTIONS:
   !
   INTEGER FUNCTION binom(n,k)

      IMPLICIT NONE

      INTEGER,       INTENT(IN)  :: n
      INTEGER,       INTENT(IN)  :: k

      binom = fac(n)/(fac(k)*fac(n-k))

   END FUNCTION binom


   ELEMENTAL REAL FUNCTION  fac(n)

      IMPLICIT NONE
 
      INTEGER, INTENT (IN) :: n
      INTEGER :: i

      fac = 0
      IF (n.LT.0) RETURN
      fac = 1
      IF (n.EQ.0) RETURN
      DO i = 2,n
      fac = fac * i
      ENDDO
 
   END FUNCTION  fac


END MODULE m_fock_basis
