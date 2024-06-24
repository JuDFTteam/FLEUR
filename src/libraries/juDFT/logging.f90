module m_juDFT_logging
    implicit none
    PRIVATE
    integer:: fh=-1

    integer,parameter:: logmode_status=1
    integer,parameter:: logmode_info=2
    integer,parameter:: logmode_warning=3
    integer,parameter:: logmode_error=4
    integer,parameter:: logmode_bug=5

    type t_log_message
        character(:),ALLOCATABLE       :: key
        character(:),allocatable       :: message
        type(t_log_message),Pointer :: next=>null()
        contains    
        procedure :: delete_all
        procedure :: add
        procedure :: write
        procedure :: report
    end type
    
    public:: t_log_message,log_start,log_stop
    public:: logmode_bug,logmode_error,logmode_warning,logmode_info,logmode_status
    contains

    recursive subroutine delete_all(logmessage)
        CLASS(t_log_message),INTENT(INOUT)::logmessage

        if (ASSOCIATED(logmessage%next)) then 
            call logmessage%next%delete_all()
            deallocate(logmessage%next)
            logmessage%next=>null()
        endif    
        if (allocated(logmessage%key)) deallocate(logmessage%key)
        if (allocated(logmessage%message)) deallocate(logmessage%message)
    end subroutine

    recursive subroutine add(logmessage,key,message)
        CLASS(t_log_message),INTENT(INOUT)::logmessage
        character(len=*),intent(in):: key,message

        if (allocated(logmessage%key)) THEN
            if (.not.associated(logmessage%next)) THEN 
                allocate(logmessage%next)
                logmessage%next%next=>null()    
            endif 
            call logmessage%next%add(key,message)
        else 
            logmessage%key=key
            logmessage%message=message
        endif 
    end subroutine     
    
    recursive subroutine write(logmessage)
        CLASS(t_log_message),INTENT(IN)::logmessage
        
        if (allocated(logmessage%key)) THEN
            if (ASSOCIATED(logmessage%next)) THEN
                write(fh,"(5a)") '  "',logmessage%key,'":"',logmessage%message,'",'
                call logmessage%next%write()
            else
                write(fh,"(5a)") '  "',logmessage%key,'":"',logmessage%message,'"'
            endif 
        endif    
    end subroutine    

    
    subroutine report(logmessage,level)
        integer,intent(in)::level
        logical::first=.true.
        CLASS(t_log_message),intent(inout):: logmessage

        integer:: dt(8)

        if (fh==-1) return
        if (first) THEN
            write(fh,"(a)") "[{"
            first=.false.
        else 
            write(fh,"(a)") ",{"
        endif
        select case(level)
        case(logmode_status)
            write(fh,"(a)") '  "loglevel":"status",'
        case(logmode_info)
            write(fh,"(a)") '  "loglevel":"info",'
        case(logmode_warning)
            write(fh,"(a)") '  "loglevel":"warning",'
        case(logmode_error)
            write(fh,"(a)") '  "loglevel":"error",'
        case(logmode_bug)
            write(fh,"(a)") '  "loglevel":"bug",'
        case default
            write(fh,"(a)") '  "loglevel":"unkown",'
        end select

        call date_and_time(values=dt)
        write(fh,"(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2,a,i2.2,a)") '  "timestamp":"',dt(3),".",dt(2),".",dt(1)," ",dt(5),":",dt(6),":",dt(7),'",'

        call logmessage%write()
        call logmessage%delete_all()

        write(fh,"(a)") "}"
        call flush(fh)
    end subroutine

    !two subroutines to start/stop logging
    subroutine log_start(filename,name)
        character(len=*),intent(in),optional:: filename,name
        type(t_log_message)::log
        fh=9231
        if (present(filename)) THEN
            open(fh,file=filename,status="replace",action="write")
        else 
            open(fh,file="runlog.json",status="replace",action="write")
        endif  
        if (present(name)) then 
            call log%add("Started",name)
        else
            call log%add("Started","FLEUR")
        endif 
        call log%report(logmode_status)    
    end subroutine 



    subroutine log_stop()
        if (fh/=-1) THEN
            write(fh,*) "]"
            close(fh)
        endif
        fh=-1    
    end subroutine
end module    
