program test_prec
    integer,PARAMETER:: d=selected_real_kind(14)
        
    if(kind(1.0)==kind(1.0_d)) stop 

    stop 1
end