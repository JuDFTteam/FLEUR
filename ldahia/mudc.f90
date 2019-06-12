MODULE m_mudc 

   CONTAINS

   SUBROUTINE mudc(U,J,l,n,l_amf,mu,jspins)

      REAL,    INTENT(IN)  :: U
      REAL,    INTENT(IN)  :: J 
      INTEGER, INTENT(IN)  :: l
      REAL,    INTENT(IN)  :: n(jspins)
      LOGICAL, INTENT(IN)  :: l_amf
      REAL,    INTENT(OUT) :: mu
      INTEGER, INTENT(IN)  :: jspins

      REAL vdcfll1,vdcfll2
      REAL vdcamf1,vdcamf2
      REAL nup,ndn

      IF(jspins.EQ.2) THEN
         nup = n(1)
         ndn = n(2)
      ELSE
         nup = 0.5 * n(1)
         ndn = nup
      ENDIF

      vdcfll1= u*(nup+ndn -0.5) - j*(nup-0.5)
      vdcfll2= u*(nup+ndn -0.5) - j*(ndn-0.5)
      vdcamf1= u*ndn+2.0*l/(2.0*l+1)*(u-j)*nup
      vdcamf2= u*nup+2.0*l/(2.0*l+1)*(u-j)*ndn
      WRITE(6,"(A)") 'Double counting chemical potential:'
      WRITE(6,9040) 'FLL: ','spin-up','spin-dn','(up+dn)/2','up-dn'
      WRITE(6,9050) vdcfll1,vdcfll2,(vdcfll1+vdcfll2)/2,vdcfll1-vdcfll2
      WRITE(6,9040) 'AMF: ','spin-up','spin-dn','(up+dn)/2','up-dn'
      WRITE(6,9050)  vdcamf1,vdcamf2,(vdcamf1+vdcamf2)/2,vdcamf1-vdcamf2

      IF(l_amf) THEN
         WRITE(6,"(A)") "Using the around-mean-field limit"
         mu = (vdcamf1+vdcamf2)/2
      ELSE
      WRITE(6,"(A)") "Using the fully-localized limit"
         mu = (vdcfll1+vdcfll2)/2
      ENDIF 
      WRITE(6,9060) mu

9040  FORMAT(TR3,A4,TR1,A7,TR3,A7,TR3,A9,TR3,A5)
9050  FORMAT(TR7,f8.4,TR2,f8.4,TR2,f8.4,TR4,f8.4)
9060  FORMAT(TR3,"mu = ",f7.4)
   END SUBROUTINE mudc

END MODULE m_mudc