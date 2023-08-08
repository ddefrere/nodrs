; ------------------------------------------------------------------------------
; ------------------------------------------------------------------------------
;+
; NAME:
;       MULTITHREAD
;
; PURPOSE:
;       distribute a command for execution across multiple cpus.
;
; EXPLANATION:
;       This procedure executes the command in the string CMD on multiple 
;       threads. Variables that must be predefined for the CMD
;       to execute must be defined in the calling procedure and specified 
;       in the string array IN_VARS. 
;
;       Variables that are to be returned from the threads are specified in
;       OUT_VARS. 
;
; CATEGORY:
;       general utilities
;
; CALLING SEQUENCE:
;       MULTITHREAD, CMD, IN_VARS [, OUT_VARS]
;
; INPUTS:
;	CMD:          A command encoded in a string, c.f. the "execute" command.
;                     the partitioning of the command across the CPUs is achieved
;                     by encoding the "[I:J]" string into the IN_VARS.
;
;	IN_VARS:      The variables from the calling scope level to transfer to
;                     the other IDLBridge sessions for inclusion in the processing.
;                     Any array variables that are to be partitioned across the CPUs 
;                     should have '[I:J]' appended in the command, e.g. 'Z=Y[I:J]'                     
;
;                     NOTE: only scalar or array variables of numeric or string type
;                     can be transferred, structures, pointers and object references
;                     cannot be transferred.
;
; KEYWORD PARAMETERS:
;
;       NTHREADS:     Specifies the number of threads to use (1..!CPU.HW_NCPU). 
;                     Default value is !CPU.HW_NCPU.
;
;       THREADS:      Returns pointers to the IDL_Bridge objects for each running
;                     process. This is only relevant if the /WAIT keyword is not set.
;
;       LOG:          If set, the output from each thread is created in the /tmp
;                     directory with the name "IDL_BRIDGE_MT_<THREAD_NUM>_<RND_NUM>.log.
;                     This is particularly useful for debugging purposes.
;
;       WAIT:         If set, then the procedure hangs, waiting for completion
;                     of all running threads.
;
;       IDL_STARTUP:  If set, the IDL_STARTUP script is run in each thread as it is
;                     created.
;
;	DONT_COPY_PATH: If set, the !PATH environment variable is NOT copied from the
;                     current IDL session to the IDLBridge sessions.
;
; OPTIONAL OUTPUTS:
;
;	OUT_VARS:     The variables from the calling scope level to be created as
;                     an output of the command.
;
;                     NOTE: only scalar or array variables of numeric or string type
;                     can be transferred, structures, pointers and object references
;                     cannot be transferred.
;
; COMMON BLOCKS:
;      Multithread2:  contains the variable with pointers to the various
;                     IDL_BRIDGE objects created.
;
; RESTRICTIONS:
;      Only works with multiple CPUs on the same machine.
;      The input/output variables can only be numeric or string scalars or arrays.
;      It is not possible to transfer structures, pointers and object references
;      between processes.
;
; PROCEDURE:
;      IDL_Bridge objects are created and the variables specified in the IN_VARS
;      string array are copied to each of the created object. Then the CMD is
;      sent to each of the Bridge objects for execution. If the /WAIT command is
;      specified, then the process waits until all objects have finished running.
;      Then the data is harvested from the processes and the output variables
;      created in the calling scope. 
;
; EXAMPLES:
;       Calculate Z = X^2 + Y^2 (a trivial example):
;          X = FINDGEN(10000L)
;          Y = X
;          multithread, 'Z=x[i:j]^2 + y[i:j]^2', ['x','y'], ['z'],/WAIT
;
;       Process multiple files:
;
;          Filenames = FILE_SEARCH('*.txt')
;          multithread, 'process_files, filenames[i:j]', ['filenames'],/WAIT
;
; MODIFICATION HISTORY:
;	Written by: H. D. R. Evans, ESA/ESTEC, Jan. 2009.
;
;
; ------------------------------------------------------------------------------
; ------------------------------------------------------------------------------
; PRO MULTITHREAD_CALLBACK - The routine called when the IDL_Bridge object finishes
;                        its task
;                                     

pro multithread_callback, status, error, objref, userdata
COMMON Multithread2, mt_obj

        ; status:
        ;    value      description             error keyword contents
        ;       0       idle                    <NA>
        ;       1       executing               <NA>
        ;       2       completed command       NULL string
        ;       3       Error halted execution  Descriptive error string
        ;       4       Aborted executeion      abort message
        ;
        ; NOTE: as this function is only called when a command ends, the
        ;       status =0 or 1 will never be returned.
        
        ;print,'Status: ',status,', id: ',userdata.threadnum, ' "',error,'"'
        ;print,'i=',userdata.start_i,', j=',userdata.stop_i

        out_vars = userdata.out_vars
        scope = userdata.scope

return
end

; ------------------------------------------------------------------------------
; ------------------------------------------------------------------------------
; PRO MULTITHREAD_WAIT - waits until the threads finish processing,
;                        then returns to the caller.
;                        DELAY=DELAY - time to wait between checks of
;                                      the status of the threads.
;                                     
pro multithread_wait,delay=delay
       COMMON Multithread2, mt_obj ; the object shared between here and the callback procedure

                                ; Now cycle, checking the status on a regular basis.

       if N_ELEMENTS(mt_obj) eq 0 then begin
           message,'No threads created, exiting',/INFO
           return
       endif

       if n_elements(delay) eq 0 then delay = 10 ; seconds
       nthread  = n_elements(mt_obj)

       status_str =['idle','running','completed','error','aborted']
       executing = intarr( nthread)
       thread_state = intarr( nthread) * 0

       still_running = 0L eq 0L

       while ( still_running ) do begin
           still_running = 0L eq 1L
                                ; loop over the threads, checking that they're still running.
           for ii=0L, nthread-1 do begin
               s = mt_obj[ii]->Status()
                                ; save the status if the previous status was "Running"
               if thread_state[ii] LE 1L then thread_state[ii]=s
               still_running = still_running OR (thread_state[ii] EQ 1L)
           endfor           
           if still_running then wait,delay
       endwhile
       
       return
   end

; ------------------------------------------------------------------------------
; ------------------------------------------------------------------------------
;
; This procedure loops over the IDL_Bridge objects, harvesting the data calculated
; from them and re-assembling the parts back into a single variable.
;
; NOTE: This doesn't work so well for multi-dimensional arrays...To now it only 
;       handles 1, 2, and 3 dimensional arrays...perhaps someone could recommend a
;       better way to do this reconstruction...

function multithread_getoutput, varname
       COMMON Multithread2, mt_obj ; the object shared between here and the callback procedure

       if n_elements(mt_obj) le 0 then begin
           message,'No threads active, returning'
           return,!values.f_nan
       endif

       mt_obj[0]->GetProperty,USERDATA=ud
       npts = ud.npts
                                ; get the data type from the first
                                ; thread, and then build variable to
                                ; hold all of the results

       d = mt_obj[0]->GetVar(varname)
       s = size( d)
       ndim = size( d,/N_DIM)
       s[s[0]] = npts
       data = MAKE_ARRAY( size=s)

       ; Now loop over the threads, harvesting the calculated values.
       for ii=0L, n_elements(mt_obj)-1 do begin 
              mt_obj[ii]->GetProperty,USERDATA=ud
              i = ud.start_i
              j = ud.stop_i
              print,'Processing ',varname, ' data from thread: ', STRTRIM(ii,2),' ',STRTRIM(i,2),'->',STRTRIM(j,2)
              if ndim le 1 then data[i:j]     = mt_obj[ii]->GetVar(varname)
              if ndim eq 2 then data[*,i:j]   = mt_obj[ii]->GetVar(varname)
              if ndim eq 3 then data[*,*,i:j] = mt_obj[ii]->GetVar(varname)
       endfor

       return,data
end

; ------------------------------------------------------------------------------
; ------------------------------------------------------------------------------
; PRO MULTITHREAD - this procedure executes the command in the string
                                ; CMD on multiple threads. Variables
                                ; that must be predefined for the CMD
                                ; to execute must exist in the calling
                                ; procedure and listed in the string
                                ; array IN_VARS. 

                                ; Variables that are to be returned
                                ; from the threads are specified in
                                ; OUT_VARS. 

                                ; THREADS returns an array of thread
                                ; objects, from which the variables
                                ; can be harvested.


pro multithread, cmd, in_vars, out_vars, nthreads=nthreads, threads=threads, log=log, $
                 wait=wait, idl_startup=idlstartup, dont_copy_path=dont_copy_path, $
                 single_thread_cmd=single_thread_cmd, force_single=force_single
COMMON Multithread2, mt_obj   ; the object shared between here and the callback procedure


;       The next bit of code handles the case where we want to run in single threaded mode,
;       e.g. if we're already running in a multithreaded process, or debugging another routine
;       (setting up and running the multiple processes can be time consuming)
        DEFSYSV,'!MULTITHREAD_PROC', exists=already_in_MT_subprocess
        if ( already_in_MT_subprocess or keyword_set(force_single)) then begin
             if N_ELEMENTS( single_thread_cmd) gt 0 THEN BEGIN
                IF (SIZE(single_thread_cmd,/TYPE) eq 7) AND $
                   STRLEN(single_thread_cmd) gt 0 $
                THEN BEGIN

                   ; We're going to try and run the command specified in the single_thread_cmd
                   ; here...

                   in_vars2432423=in_vars           ; we suffix the main variables to prevent
                   cmd2432423 = single_thread_cmd   ; name space collisions....
                   out_vars2432423 = out_vars

                   ; we should also check that there aren't any in_vars that clash with the
                   ; already unlikely variable names in this section of code. TBD

                   ; Copy in the required variables from the calling scope
                   FOR iiii2432423=0L,n_elements(in_vars2432423)-1 do begin
                        ; copy in the input variables to this scope
                        (scope_varfetch( in_vars2432423[iiii2432423], level=0,/ENTER)) = $
                            SCOPE_VARFETCH( in_vars2432423[iiii2432423], LEVEL=-1,/ENTER)
                   ENDFOR
                   ; Run the command
                   result2432423 = EXECUTE( cmd2432423 )
                   ; copy the resulting variables into the calling scope
                   FOR i2432423=0L,n_elements(out_vars2432423)-1 do begin
                       if strlen(out_vars2432423[i2432423]) gt 0 then begin
                           (scope_varfetch(out_vars2432423[i2432423],LEVEL=-1,/ENTER))= $
                                SCOPE_VARFETCH(out_vars2432423[i2432423],LEVEL=0,/ENTER)
                       ENDIF
                   ENDFOR
                ENDIF ELSE $
                    message,'Single threaded command is not a string type.'
             ENDIF ELSE $
                message,'Already in a multithreaded subprocess and no single_thread_cmd specified.'
             RETURN
        endif

        curr_scope = scope_level()
        prev_scope = curr_scope - 1

        avail_vars = scope_varname(level=prev_scope)
        
                                ; We assume that the number of points
                                ; to calculate is that specified by
                                ; the variable named by the IN_VARS[0]
                                ; parameter.
        npts = n_elements( scope_varfetch(in_vars[0],level=prev_scope))

        if n_elements(nthreads) eq 0 then nthreads=!cpu.hw_ncpu
        nthreads = (nthreads > 1) < MIN([!CPU.HW_NCPU, npts])

        mt_obj = objarr(nthreads)

        if nthreads gt 1 then begin
                start_i = LINDGEN(nthreads) * npts / nthreads
                stop_i  = [start_i[1:*] - 1L, npts-1]
        endif else begin
                start_i = 0L
                stop_i  = npts-1
        endelse

        for i=0,nthreads-1 do begin

                                ; Create the new thread if it doesn't exist
            if NOT (PTR_VALID( mt_obj[i])) THEN begin
                message,/INFO,'Creating process: '+STRTRIM(i,2)
                if keyword_set(log) then begin
                    output = '/tmp/IDL_Bridge_MT_' + $
                              STRING(i,form='(i3.3)') + '_' + $
                              STRTRIM(STRING(SYSTIME(/SEC)*1000.0,FORM='(Z)'),2)
                    mt_obj[i] = OBJ_NEW('IDL_IDLBridge',$
                                    callback='multithread_callback',$
                                    output=output+'.log')
                endif else begin
                    mt_obj[i] = OBJ_NEW('IDL_IDLBridge',$
                                    callback='multithread_callback')
                endelse
            endif

            ; run the IDL_STARTUP command to initialise the thread.
            if keyword_set( idl_startup) then $
               mt_obj[i]->EXECUTE,'@' + PREF_GET('IDL_STARTUP')

                                ; setup the state and user variables
            if n_elements(out_vars) gt 0 then $
               state = {oBridge:mt_obj[i], $
                     threadNum:i, $
                     start_i:0L, stop_i:0L, NPTS:0L, $
                     scope:prev_scope, $
                     out_vars:out_vars} $
            else $
               state = {oBridge:mt_obj[i], $
                     threadNum:i, $
                     start_i:0L, stop_i:0L, NPTS:0L, $
                     scope:prev_scope, $
                     out_vars:0L}

            state.start_i = start_i[i]   ; Set the index at which to start processing
            state.stop_i = stop_i[i]     ; Set the index at which to stop processing
            state.NPTS   = npts          ; Save the total number of points to process

            mt_obj[i]->setproperty,userdata=state

            ; Set up a system variable in the slave process to allow us to 
            ; detect when a program is running in a slave process. This can
            ; prevent the nested running of this system, e.g. a slave process
            ; itself tries to use this system to multithread. We only want
            ; one level of multithreading.
            mt_obj[i]->Execute,"DEFSYSV,'!MULTITHREAD_PROC', exists=in_mt"
            mt_obj[i]->Execute,"IF ( ~in_mt) then DEFSYSV,'!MULTITHREAD_PROC', 1, 1"

                                ; load up shared variables in the new thread
            for j=0,n_elements(in_vars)-1 do begin
                                ; copy the variables from the calling
                                ; procedure to the new thread's memory
                                ; space
                mt_obj[i]->SetVar,in_vars[j],(scope_varfetch(in_vars[j],LEVEL=-1))
            endfor
                                ; copy/set the index variables from
                                ; that are set in the CMD function,
                                ; and ensure the path is set
            mt_obj[i]->SetVar,'i',start_i[i]
            mt_obj[i]->SetVar,'j',stop_i[i]
            cd,current=cwd
            if NOT KEYWORD_SET( DONT_COPY_PATH) then begin
               this_path = cwd + PATH_SEP(/SEARCH) + !PATH
               mt_obj[i]->SetVar,'this_path',this_path
               mt_obj[i]->Execute,'!PATH=this_path'
            endif
                                ; Execute the command on the thread,
                                ; but don't wait around for completion.
            mt_obj[i]->Execute,'cd,"'+cwd+'"'
            mt_obj[i]->Execute,cmd,/NOWAIT
        endfor
        

        ; wait here until all the threads finish...
        if keyword_set(wait) then BEGIN
                                ; wait until finished
            multithread_wait

                                ; loop over output variables, creating
                                ; them in the calling scope and
                                ; harvesting the results from the threads.

            IF N_ELEMENTS(out_vars) gt 0 THEN BEGIN
               FOR i=0L,n_elements(out_vars)-1 do begin
                   if strlen(out_vars[i]) gt 0 then $
                   (scope_varfetch(out_vars[i],LEVEL=-1,/ENTER)) = multithread_getoutput(out_vars[i])
               ENDFOR
            ENDIF

            ; clean up the threads, they've served their purpose.

;            for i=0L,n_elements(mt_obj)-1 do begin
                ;print,'Killing thread #',i
                ;obj_destroy,mt_obj[i]
;            endfor

        ENDIF ELSE threads= mt_obj  ; return the threads to the calling routine
        return
    end
    
