	  9  7   k820309              12.1        f��[                                                                                                           
       mcio.f90 MCIO                                                     
       LINELEN CHKOK CHKERR CHKEOF                                                                                                    256                                                                                                    0                                                                                                   1                                                                                       ��������        %         @                                                           #MCIO_GETUNIT%PRESENT    #U_TRY                                                    PRESENT           
 @                                          $        @                                	                           #MCIO_GETCAPITAL%ACHAR 
   #MCIO_GETCAPITAL%IACHAR    #MCIO_GETCAPITAL%LEN_TRIM    #MCIO_GETCAPITAL%LEN    #CHR_IN    H r      5 O p                                                               
     ACHAR                                                 IACHAR                                                 LEN_TRIM                                                 LEN           
  @                                                  1 $        @                                                           #MCIO_CENTER%REPEAT    #MCIO_CENTER%TRIM    #MCIO_CENTER%ADJUSTL    #MCIO_CENTER%LEN_TRIM    #CHR    #N    5 O p                                                        REPEAT                                                 TRIM                                                 ADJUSTL                                                 LEN_TRIM           
  @                                                  1           
                                             %         @                                                           #MCIO_ISODD%MOD    #VAL                                                    MOD           
  @                                          %         @                                                           #MCIO_ASSERT_FILE%TRIM    #MCIO_ASSERT_FILE%ADJUSTL    #FNAME                                                    TRIM                                                 ADJUSTL           
  @                                                  1 #         @                                                    #MCIO_FINDLINE%TRIM    #MCIO_FINDLINE%ADJUSTL    #U_ID     #CHR_IN !   #CHR_OUT "   #IERR #                                                   TRIM                                                 ADJUSTL           
                                                        
  @                              !                    1           D                                "                     1           D                                 #            #         @                                  $                  #MCIO_PRINTHEAD%REPEAT %   #MCIO_PRINTHEAD%PRESENT &   #CHR_IN '   #SYMB (   #U )                                              %     REPEAT                                            &     PRESENT           
  @                              '                    1           
 @                               (                                     
 @                               )           #         @                                  *                  #MCIO_PRINTERR%TRIM +   #MCIO_PRINTERR%ADJUSTL ,   #MCIO_PRINTERR%PRESENT -   #ACT .   #CHR_IN /   #LINE 0                                              +     TRIM                                            ,     ADJUSTL                                            -     PRESENT           
  @                              .                    1           
  @                              /                    1           
 @                              0                    1 %         @                                 1                          #MCIO_FINDSTR%TRIM 2   #MCIO_FINDSTR%ADJUSTL 3   #MCIO_FINDSTR%LEN 4   #LINE 5   #STR 6                                              2     TRIM                                            3     ADJUSTL                                            4     LEN           
  @                              5                    1           
  @                              6                    1    �         fn#fn    �   \   J  CONSTANTS "     s       LINELEN+CONSTANTS     �  q       CHKOK+CONSTANTS !   �  q       CHKERR+CONSTANTS !   g  p       CHKEOF+CONSTANTS    �  u       MCIO_GETUNIT %   L  @      MCIO_GETUNIT%PRESENT #   �  @   a   MCIO_GETUNIT%U_TRY     �  �       MCIO_GETCAPITAL &   �  >      MCIO_GETCAPITAL%ACHAR '     ?      MCIO_GETCAPITAL%IACHAR )   G  A      MCIO_GETCAPITAL%LEN_TRIM $   �  <      MCIO_GETCAPITAL%LEN '   �  L   a   MCIO_GETCAPITAL%CHR_IN      �       MCIO_CENTER #   �  ?      MCIO_CENTER%REPEAT !   $  =      MCIO_CENTER%TRIM $   a  @      MCIO_CENTER%ADJUSTL %   �  A      MCIO_CENTER%LEN_TRIM     �  L   a   MCIO_CENTER%CHR    .  @   a   MCIO_CENTER%N    n  m       MCIO_ISODD    �  <      MCIO_ISODD%MOD    	  @   a   MCIO_ISODD%VAL !   W	  �       MCIO_ASSERT_FILE &   �	  =      MCIO_ASSERT_FILE%TRIM )   (
  @      MCIO_ASSERT_FILE%ADJUSTL '   h
  L   a   MCIO_ASSERT_FILE%FNAME    �
  �       MCIO_FINDLINE #   \  =      MCIO_FINDLINE%TRIM &   �  @      MCIO_FINDLINE%ADJUSTL #   �  @   a   MCIO_FINDLINE%U_ID %     L   a   MCIO_FINDLINE%CHR_IN &   e  L   a   MCIO_FINDLINE%CHR_OUT #   �  @   a   MCIO_FINDLINE%IERR    �  �       MCIO_PRINTHEAD &   �  ?      MCIO_PRINTHEAD%REPEAT '   �  @      MCIO_PRINTHEAD%PRESENT &     L   a   MCIO_PRINTHEAD%CHR_IN $   X  P   a   MCIO_PRINTHEAD%SYMB !   �  @   a   MCIO_PRINTHEAD%U    �  �       MCIO_PRINTERR #   �  =      MCIO_PRINTERR%TRIM &   �  @      MCIO_PRINTERR%ADJUSTL &     @      MCIO_PRINTERR%PRESENT "   Z  L   a   MCIO_PRINTERR%ACT %   �  L   a   MCIO_PRINTERR%CHR_IN #   �  L   a   MCIO_PRINTERR%LINE    >  �       MCIO_FINDSTR "   �  =      MCIO_FINDSTR%TRIM %   %  @      MCIO_FINDSTR%ADJUSTL !   e  <      MCIO_FINDSTR%LEN "   �  L   a   MCIO_FINDSTR%LINE !   �  L   a   MCIO_FINDSTR%STR 