	  �Z  �   k820309              12.1        ؿ�[                                                                                                           
       sfsetup.f90 SFSETUP              READ_SETUP_PARAMETERS LOAD_SETUP LOAD_SETUP_ASCII SAVE_SETUP SAVE_SETUP_ASCII SKIP_SETUP NEW_SETUP DEL_SETUP STP_INIT STP_FINAL STP_GET_RANGE STP_EVAL STP_NORMALIZE STP_NSF_MAX STP_PRINT_INFO STP_ASSERT_INIT STP_ASSERT_MODULEINIT SETUP SFVAL SFDERIV_I SFDERIV_J NSF_MAX NNB_MAX                      @                              
       AEIO_READLINE TYPELEN                      @                              
       IO_LOWER IO_READVAL IO_ADJUSTL IO_UNIT                      @                              
       SF_INIT SF_FINAL SF_ADD_RAD SF_ADD_ANG SF_FINGERPRINT                �                                      u #AEIO_READLINE_C1    #AEIO_READLINE_CN    #AEIO_READLINE_I1    #AEIO_READLINE_IN    #AEIO_READLINE_D1 !   #AEIO_READLINE_DN '   #AEIO_READLINE_L1 .   #AEIO_READLINE_LN 4   #         @     @                                            #AEIO_READLINE_C1%PRESENT    #AEIO_READLINE_C1%TRIM    #AEIO_READLINE_C1%ADJUSTL    #AEIO_READLINE_C1%LEN_TRIM    #U_IN 	   #ILINE 
   #LINE    #STAT                                                    PRESENT                                                 TRIM                                                 ADJUSTL                                                 LEN_TRIM           
                                  	                     
                                 
                                                                           1                                                        #         @     @                                             #AEIO_READLINE_CN%PRESENT    #U_IN    #ILINE    #LINE    #N    #STAT                                                    PRESENT           
                                                       
                                             ,                                                                   p          5 O p            5 O p                          1           
                                                                                                    #         @     @                                             #AEIO_READLINE_I1%PRESENT    #U_IN    #ILINE    #LINE    #STAT                                                    PRESENT           
                                                       
                                                                                                                                                           #         @     @                                             #AEIO_READLINE_IN%PRESENT    #U_IN    #ILINE    #LINE    #N    #STAT                                                     PRESENT           
                                                       
                                                                                                                p          5 O p            5 O p                                    
                                                                                                     #         @     @                           !                  #AEIO_READLINE_D1%PRESENT "   #U_IN #   #ILINE $   #LINE %   #STAT &                                              "     PRESENT           
                                  #                     
                                 $                                                       %     
                                                  &            #         @     @                           '                  #AEIO_READLINE_DN%PRESENT (   #U_IN )   #ILINE *   #LINE +   #N ,   #STAT -                                              (     PRESENT           
                                  )                     
                                 *                                                      +                    
     p          5 O p            5 O p                                    
                                  ,                                                      -            #         @     @                           .                  #AEIO_READLINE_L1%PRESENT /   #U_IN 0   #ILINE 1   #LINE 2   #STAT 3                                              /     PRESENT           
                                  0                     
                                 1                                                       2                                                       3            #         @     @                           4                  #AEIO_READLINE_LN%PRESENT 5   #U_IN 6   #ILINE 7   #LINE 8   #N 9   #STAT :                                              5     PRESENT           
                                  6                     
                                 7                                                      8                         p          5 O p            5 O p                                    
                                  9                                                      :                           �                                       u #ADJUSTL_I ;   #ADJUSTL_D ?   $         @    @                          ;     2                     #ADJUSTL_I%TRIM <   #ADJUSTL_I%ADJUSTL =   #INT >                         @                            <     TRIM               @                            =     ADJUSTL           
   @                              >           $         @    @                          ?     2                     #ADJUSTL_D%TRIM @   #ADJUSTL_D%ADJUSTL A   #ADJUSTL_D%PRESENT B   #DP C   #DIGITS D                         @                            @     TRIM               @                            A     ADJUSTL               @                            B     PRESENT           
   @                              C     
                
  @                              D                          �                                      u #READVAL_I1 E   #READVAL_IN L   #READVAL_D1 T   #READVAL_DN [   #READVAL_C1 c   #READVAL_L j   #         @     @                          E                  #READVAL_I1%TRIM F   #READVAL_I1%LEN_TRIM G   #READVAL_I1%INDEX H   #STRING I   #NAME J   #VAL K                 @                            F     TRIM               @                            G     LEN_TRIM               @                            H     INDEX           
   @                             I                    1           
   @                             J                    1           
  @                              K            #         @     @                           L                  #READVAL_IN%TRIM M   #READVAL_IN%LEN_TRIM N   #READVAL_IN%INDEX O   #STRING P   #NAME Q   #VAL R   #N S                 @                            M     TRIM               @                            N     LEN_TRIM               @                            O     INDEX           
   @                             P                    1           
   @                             Q                    1          
  @                              R                         p          5 O p            5 O p                                    
   @                              S           #         @     @                          T                  #READVAL_D1%TRIM U   #READVAL_D1%LEN_TRIM V   #READVAL_D1%INDEX W   #STRING X   #NAME Y   #VAL Z                 @                            U     TRIM               @                            V     LEN_TRIM               @                            W     INDEX           
   @                             X                    1           
   @                             Y                    1           
  @                              Z     
       #         @     @                           [                  #READVAL_DN%TRIM \   #READVAL_DN%LEN_TRIM ]   #READVAL_DN%INDEX ^   #STRING _   #NAME `   #VAL a   #N b                 @                            \     TRIM               @                            ]     LEN_TRIM               @                            ^     INDEX           
   @                             _                    1           
   @                             `                    1          
  @                              a                    
     p          5 O p            5 O p                                    
   @                              b           #         @     @                          c                  #READVAL_C1%TRIM d   #READVAL_C1%LEN_TRIM e   #READVAL_C1%INDEX f   #STRING g   #NAME h   #VAL i                 @                            d     TRIM               @                            e     LEN_TRIM               @                            f     INDEX           
   @                             g                    1           
   @                             h                    1           
  @                             i                     1 #         @     @                           j                  #READVAL_L%TRIM k   #READVAL_L%INDEX l   #STRING m   #NAME n   #VAL o                 @                            k     TRIM               @                            l     INDEX           
   @                             m                    1           
   @                             n                    1             @                              o                                                         p                                                      2$        @                                q                           #IO_LOWER%LEN r   #IO_LOWER%ICHAR s   #IO_LOWER%CHAR t   #STR_IN u   H r r     5 O p                                  @                            r     LEN               @                            s     ICHAR               @                            t     CHAR           
   @                             u                    1 %         @                                v                          #IO_UNIT%PRESENT w   #U_TRY x                 @                            w     PRESENT           
  @                              x           &         @                                 y     �                    #READ_SETUP_PARAMETERS%LEN_TRIM z   #READ_SETUP_PARAMETERS%INDEX {   #READ_SETUP_PARAMETERS%TRIM |   #READ_SETUP_PARAMETERS%SIZE }   #READ_SETUP_PARAMETERS%LEN ~   #INFILE    #GLOBAL_TYPES �   #SETUP �                 @                            z     LEN_TRIM               @                            {     INDEX               @                            |     TRIM               @                            }     SIZE               @                            ~     LEN           
  @@                                                 1 ,          
 @@                             �                                  &                                           1 &         @                                 �     �                    #LOAD_SETUP%PRESENT �   #LOAD_SETUP%TRIM �   #LOAD_SETUP%SIZE �   #GLOBAL_TYPES �   #FILE �   #UNIT �   #SETUP �                 @                            �     PRESENT               @                            �     TRIM               @                            �     SIZE ,          
 @@                             �                                  &                                           1           
 @@                             �                    1           
 @@                              �           &         @                                 �     �                    #LOAD_SETUP_ASCII%PRESENT �   #LOAD_SETUP_ASCII%TRIM �   #LOAD_SETUP_ASCII%SIZE �   #GLOBAL_TYPES �   #FILE �   #UNIT �   #SETUP �                 @                            �     PRESENT               @                            �     TRIM               @                            �     SIZE ,          
 @@                             �                                  &                                           1           
 @@                             �                    1           
 @@                              �           #         @                                  �                  #SAVE_SETUP%PRESENT �   #SAVE_SETUP%TRIM �   #STP �   #FILE �   #UNIT �                 @                            �     PRESENT               @                            �     TRIM           
  @@                              �     �             #SETUP �             
 @@                             �                    1           
 @@                              �           #         @                                  �                  #SAVE_SETUP_ASCII%PRESENT �   #SAVE_SETUP_ASCII%TRIM �   #STP �   #FILE �   #UNIT �                 @                            �     PRESENT               @                            �     TRIM           
  @@                              �     �             #SETUP �             
 @@                             �                    1           
 @@                              �           #         @                                  �                   #U �             
   @                              �           &         @                                 �     �                     #NSF �   #NENV �   #NTYPES_GLOBAL �   #SETUP �             
   @                              �                     
   @                              �                     
   @                              �           #         @                                  �                   #STP �             
D  @                              �     �              #SETUP �   #         @                                  �                  #STP_INIT%MAX �   #STP_INIT%TRIM �   #NTYPES �   #STP �   #N_NB_MAX �                 @                            �     MAX               @                            �     TRIM           
  @@                              �                    
  @@                              �             �           p          5 � p        r �       5 � p        r �                     #SETUP �             
   @                              �           #         @                                  �                  #STP_FINAL%TRIM �   #NTYPES �   #STP �                 @                            �     TRIM           
   @                              �                    
   @                              �             �           p          5 � p        r �       5 � p        r �                     #SETUP �   #         @                                  �                  #STP_GET_RANGE%MIN �   #STP_GET_RANGE%MAX �   #NTYPES �   #STP �   #RC_MIN �   #RC_MAX �                 @                            �     MIN               @                            �     MAX           
   @                              �                    
  @@                              �             �           p          5 � p        r �       5 � p        r �                     #SETUP �             D @@                              �     
                 D @@                              �     
       #         @                                  �                  #STP_EVAL%DBLE �   #STP_EVAL%MIN �   #STP_EVAL%MAX �   #STP_EVAL%PRESENT �   #STP_EVAL%TRIM �   #ITYPE0 �   #COO0 �   #N �   #COO1 �   #TYPE1 �   #STP �   #DERIV �   #SCALED �                 @                            �     DBLE               @                            �     MIN               @                            �     MAX               @                            �     PRESENT               @                            �     TRIM           
   @                              �                     
  @@                              �                   
    p          p            p                                    
  @@                              �                    
  @@                              �                    
    p          p          5 � p        r �       p          5 � p        r �                              
   @                              �                        p          5 � p        r �       5 � p        r �                               
D @@                              �     �              #SETUP �             
 @@                              �                     
 @@                              �           #         @                                 �                  #STP_NORMALIZE%SQRT �   #STP_NORMALIZE%ABS �   #STP_NORMALIZE%PRESENT �   #STP_NORMALIZE%ADJUSTL �   #STP_NORMALIZE%TRIM �   #STP �   #DERIV �                 @                            �     SQRT               @                            �     ABS               @                            �     PRESENT               @                            �     ADJUSTL               @                            �     TRIM           
D @@                              �     �              #SETUP �             
 @@                              �           %         @                                �                          #STP_NSF_MAX%MAX �   #STP_NSF_MAX%PRESENT �   #STP_NSF_MAX%SIZE �   #STP �                 @                            �     MAX               @                            �     PRESENT               @                            �     SIZE           
@@                              �            �                     &                                           #SETUP �   #         @                                  �                  #STP_PRINT_INFO%ADJUSTL �   #STP_PRINT_INFO%TRIM �   #STP �                 @                            �     ADJUSTL               @                            �     TRIM           
@ @@                              �     �             #SETUP �   #         @                                 �                   #STP �             
   @                              �     �             #SETUP �   #         @                                 �                                      @               @           �     '�                   #INIT �   #NEVAL �   #DESCRIPTION �   #ATOMTYPE �   #NENV �   #ENVTYPES �   #NTYPES_GLOBAL �   #GTYPE �   #LTYPE �   #RC_MIN �   #RC_MAX �   #SFTYPE �   #NSF �   #SF �   #NSFPARAM �   #SFPARAM �   #SFENV �   #SFVAL_MIN �   #SFVAL_MAX �   #SFVAL_AVG �   #SFVAL_COV �                � $                              �                                � $                              �                               � $                             �                                       � $                             �                                      � $                              �                 .           � $                             �                                        &                                                                � $                              �     X                       � $                              �            `                            &                                                      � $                              �            �             	               &                                                        � $                              �     �      
   
                � $                              �     �         
                � $                             �     d                                  � $                              �     d                       � $                              �            h                            &                                                        � $                              �     �                       � $                              �            �                
            &                   &                                                      � $                              �                                        &                   &                                                      � $                              �            x                
            &                                                      � $                              �            �                
            &                                                      � $                              �                            
            &                                                      � $                              �            P                
            &                                                    @@                               �                   
                &                                                    @@                               �                   
                &                   &                                                    @@                               �                   
                &                   &                   &                                                      @                                �                       @                                �               �         fn#fn    �   &  b   uapp(SFSETUP    �  V   J  AEIO    8  g   J  IO    �  v   J  SYMMFUNC '     �       gen@AEIO_READLINE+AEIO &     �      AEIO_READLINE_C1+AEIO 6   �  @      AEIO_READLINE_C1%PRESENT+AEIO=PRESENT 0   ,  =      AEIO_READLINE_C1%TRIM+AEIO=TRIM 6   i  @      AEIO_READLINE_C1%ADJUSTL+AEIO=ADJUSTL 8   �  A      AEIO_READLINE_C1%LEN_TRIM+AEIO=LEN_TRIM +   �  @   a   AEIO_READLINE_C1%U_IN+AEIO ,   *  @   a   AEIO_READLINE_C1%ILINE+AEIO +   j  L   a   AEIO_READLINE_C1%LINE+AEIO +   �  @   a   AEIO_READLINE_C1%STAT+AEIO &   �  �      AEIO_READLINE_CN+AEIO 6   �  @      AEIO_READLINE_CN%PRESENT+AEIO=PRESENT +   �  @   a   AEIO_READLINE_CN%U_IN+AEIO ,     @   a   AEIO_READLINE_CN%ILINE+AEIO +   L  �   a   AEIO_READLINE_CN%LINE+AEIO (   �  @   a   AEIO_READLINE_CN%N+AEIO +   4	  @   a   AEIO_READLINE_CN%STAT+AEIO &   t	  �      AEIO_READLINE_I1+AEIO 6   
  @      AEIO_READLINE_I1%PRESENT+AEIO=PRESENT +   C
  @   a   AEIO_READLINE_I1%U_IN+AEIO ,   �
  @   a   AEIO_READLINE_I1%ILINE+AEIO +   �
  @   a   AEIO_READLINE_I1%LINE+AEIO +     @   a   AEIO_READLINE_I1%STAT+AEIO &   C  �      AEIO_READLINE_IN+AEIO 6   �  @      AEIO_READLINE_IN%PRESENT+AEIO=PRESENT +     @   a   AEIO_READLINE_IN%U_IN+AEIO ,   Y  @   a   AEIO_READLINE_IN%ILINE+AEIO +   �  �   a   AEIO_READLINE_IN%LINE+AEIO (   =  @   a   AEIO_READLINE_IN%N+AEIO +   }  @   a   AEIO_READLINE_IN%STAT+AEIO &   �  �      AEIO_READLINE_D1+AEIO 6   L  @      AEIO_READLINE_D1%PRESENT+AEIO=PRESENT +   �  @   a   AEIO_READLINE_D1%U_IN+AEIO ,   �  @   a   AEIO_READLINE_D1%ILINE+AEIO +     @   a   AEIO_READLINE_D1%LINE+AEIO +   L  @   a   AEIO_READLINE_D1%STAT+AEIO &   �  �      AEIO_READLINE_DN+AEIO 6   "  @      AEIO_READLINE_DN%PRESENT+AEIO=PRESENT +   b  @   a   AEIO_READLINE_DN%U_IN+AEIO ,   �  @   a   AEIO_READLINE_DN%ILINE+AEIO +   �  �   a   AEIO_READLINE_DN%LINE+AEIO (   �  @   a   AEIO_READLINE_DN%N+AEIO +   �  @   a   AEIO_READLINE_DN%STAT+AEIO &     �      AEIO_READLINE_L1+AEIO 6   �  @      AEIO_READLINE_L1%PRESENT+AEIO=PRESENT +   �  @   a   AEIO_READLINE_L1%U_IN+AEIO ,     @   a   AEIO_READLINE_L1%ILINE+AEIO +   U  @   a   AEIO_READLINE_L1%LINE+AEIO +   �  @   a   AEIO_READLINE_L1%STAT+AEIO &   �  �      AEIO_READLINE_LN+AEIO 6   k  @      AEIO_READLINE_LN%PRESENT+AEIO=PRESENT +   �  @   a   AEIO_READLINE_LN%U_IN+AEIO ,   �  @   a   AEIO_READLINE_LN%ILINE+AEIO +   +  �   a   AEIO_READLINE_LN%LINE+AEIO (   �  @   a   AEIO_READLINE_LN%N+AEIO +     @   a   AEIO_READLINE_LN%STAT+AEIO "   O  ^       gen@IO_ADJUSTL+IO    �  �      ADJUSTL_I+IO '   9  =      ADJUSTL_I%TRIM+IO=TRIM -   v  @      ADJUSTL_I%ADJUSTL+IO=ADJUSTL !   �  @   e   ADJUSTL_I%INT+IO    �  �      ADJUSTL_D+IO '   �  =      ADJUSTL_D%TRIM+IO=TRIM -   �  @      ADJUSTL_D%ADJUSTL+IO=ADJUSTL -   !  @      ADJUSTL_D%PRESENT+IO=PRESENT     a  @   e   ADJUSTL_D%DP+IO $   �  @   e   ADJUSTL_D%DIGITS+IO "   �  �       gen@IO_READVAL+IO    �  �      READVAL_I1+IO (   +  =      READVAL_I1%TRIM+IO=TRIM 0   h  A      READVAL_I1%LEN_TRIM+IO=LEN_TRIM *   �  >      READVAL_I1%INDEX+IO=INDEX %   �  L   e   READVAL_I1%STRING+IO #   3  L   e   READVAL_I1%NAME+IO "     @   e   READVAL_I1%VAL+IO    �  �      READVAL_IN+IO (   q  =      READVAL_IN%TRIM+IO=TRIM 0   �  A      READVAL_IN%LEN_TRIM+IO=LEN_TRIM *   �  >      READVAL_IN%INDEX+IO=INDEX %   -  L   e   READVAL_IN%STRING+IO #   y  L   e   READVAL_IN%NAME+IO "   �  �   e   READVAL_IN%VAL+IO     i  @   e   READVAL_IN%N+IO    �  �      READVAL_D1+IO (   T   =      READVAL_D1%TRIM+IO=TRIM 0   �   A      READVAL_D1%LEN_TRIM+IO=LEN_TRIM *   �   >      READVAL_D1%INDEX+IO=INDEX %   !  L   e   READVAL_D1%STRING+IO #   \!  L   e   READVAL_D1%NAME+IO "   �!  @   e   READVAL_D1%VAL+IO    �!  �      READVAL_DN+IO (   �"  =      READVAL_DN%TRIM+IO=TRIM 0   �"  A      READVAL_DN%LEN_TRIM+IO=LEN_TRIM *   #  >      READVAL_DN%INDEX+IO=INDEX %   V#  L   e   READVAL_DN%STRING+IO #   �#  L   e   READVAL_DN%NAME+IO "   �#  �   e   READVAL_DN%VAL+IO     �$  @   e   READVAL_DN%N+IO    �$  �      READVAL_C1+IO (   }%  =      READVAL_C1%TRIM+IO=TRIM 0   �%  A      READVAL_C1%LEN_TRIM+IO=LEN_TRIM *   �%  >      READVAL_C1%INDEX+IO=INDEX %   9&  L   e   READVAL_C1%STRING+IO #   �&  L   e   READVAL_C1%NAME+IO "   �&  L   e   READVAL_C1%VAL+IO    '  �      READVAL_L+IO '   �'  =      READVAL_L%TRIM+IO=TRIM )   �'  >      READVAL_L%INDEX+IO=INDEX $   ((  L   e   READVAL_L%STRING+IO "   t(  L   e   READVAL_L%NAME+IO !   �(  @   e   READVAL_L%VAL+IO     )  q       TYPELEN+AEIO    q)  �       IO_LOWER+IO $   :*  <      IO_LOWER%LEN+IO=LEN (   v*  >      IO_LOWER%ICHAR+IO=ICHAR &   �*  =      IO_LOWER%CHAR+IO=CHAR #   �*  L   e   IO_LOWER%STR_IN+IO    =+  p       IO_UNIT+IO +   �+  @      IO_UNIT%PRESENT+IO=PRESENT !   �+  @   e   IO_UNIT%U_TRY+IO &   -,        READ_SETUP_PARAMETERS /   J-  A      READ_SETUP_PARAMETERS%LEN_TRIM ,   �-  >      READ_SETUP_PARAMETERS%INDEX +   �-  =      READ_SETUP_PARAMETERS%TRIM +   .  =      READ_SETUP_PARAMETERS%SIZE *   C.  <      READ_SETUP_PARAMETERS%LEN -   .  L   a   READ_SETUP_PARAMETERS%INFILE 3   �.  �   a   READ_SETUP_PARAMETERS%GLOBAL_TYPES    [/  �       LOAD_SETUP #   0  @      LOAD_SETUP%PRESENT     ^0  =      LOAD_SETUP%TRIM     �0  =      LOAD_SETUP%SIZE (   �0  �   a   LOAD_SETUP%GLOBAL_TYPES     h1  L   a   LOAD_SETUP%FILE     �1  @   a   LOAD_SETUP%UNIT !   �1  �       LOAD_SETUP_ASCII )   �2  @      LOAD_SETUP_ASCII%PRESENT &   	3  =      LOAD_SETUP_ASCII%TRIM &   F3  =      LOAD_SETUP_ASCII%SIZE .   �3  �   a   LOAD_SETUP_ASCII%GLOBAL_TYPES &   4  L   a   LOAD_SETUP_ASCII%FILE &   _4  @   a   LOAD_SETUP_ASCII%UNIT    �4  �       SAVE_SETUP #   15  @      SAVE_SETUP%PRESENT     q5  =      SAVE_SETUP%TRIM    �5  S   a   SAVE_SETUP%STP     6  L   a   SAVE_SETUP%FILE     M6  @   a   SAVE_SETUP%UNIT !   �6  �       SAVE_SETUP_ASCII )   +7  @      SAVE_SETUP_ASCII%PRESENT &   k7  =      SAVE_SETUP_ASCII%TRIM %   �7  S   a   SAVE_SETUP_ASCII%STP &   �7  L   a   SAVE_SETUP_ASCII%FILE &   G8  @   a   SAVE_SETUP_ASCII%UNIT    �8  O       SKIP_SETUP    �8  @   a   SKIP_SETUP%U    9  �       NEW_SETUP    �9  @   a   NEW_SETUP%NSF    �9  @   a   NEW_SETUP%NENV (   :  @   a   NEW_SETUP%NTYPES_GLOBAL    W:  Q       DEL_SETUP    �:  S   a   DEL_SETUP%STP    �:  �       STP_INIT    �;  <      STP_INIT%MAX    �;  =      STP_INIT%TRIM     <  @   a   STP_INIT%NTYPES    D<  �   a   STP_INIT%STP "   =  @   a   STP_INIT%N_NB_MAX    C=  q       STP_FINAL    �=  =      STP_FINAL%TRIM !   �=  @   a   STP_FINAL%NTYPES    1>  �   a   STP_FINAL%STP    �>  �       STP_GET_RANGE "   �?  <      STP_GET_RANGE%MIN "   �?  <      STP_GET_RANGE%MAX %   @  @   a   STP_GET_RANGE%NTYPES "   K@  �   a   STP_GET_RANGE%STP %   
A  @   a   STP_GET_RANGE%RC_MIN %   JA  @   a   STP_GET_RANGE%RC_MAX    �A  �       STP_EVAL    �B  =      STP_EVAL%DBLE    �B  <      STP_EVAL%MIN    �B  <      STP_EVAL%MAX !   9C  @      STP_EVAL%PRESENT    yC  =      STP_EVAL%TRIM     �C  @   a   STP_EVAL%ITYPE0    �C  �   a   STP_EVAL%COO0    �D  @   a   STP_EVAL%N    �D  �   a   STP_EVAL%COO1    �E  �   a   STP_EVAL%TYPE1    RF  S   a   STP_EVAL%STP    �F  @   a   STP_EVAL%DERIV     �F  @   a   STP_EVAL%SCALED    %G  �       STP_NORMALIZE #   �G  =      STP_NORMALIZE%SQRT "   ;H  <      STP_NORMALIZE%ABS &   wH  @      STP_NORMALIZE%PRESENT &   �H  @      STP_NORMALIZE%ADJUSTL #   �H  =      STP_NORMALIZE%TRIM "   4I  S   a   STP_NORMALIZE%STP $   �I  @   a   STP_NORMALIZE%DERIV    �I  �       STP_NSF_MAX     dJ  <      STP_NSF_MAX%MAX $   �J  @      STP_NSF_MAX%PRESENT !   �J  =      STP_NSF_MAX%SIZE     K  �   a   STP_NSF_MAX%STP    �K  �       STP_PRINT_INFO '   :L  @      STP_PRINT_INFO%ADJUSTL $   zL  =      STP_PRINT_INFO%TRIM #   �L  S   a   STP_PRINT_INFO%STP     
M  Q       STP_ASSERT_INIT $   [M  S   a   STP_ASSERT_INIT%STP &   �M  H       STP_ASSERT_MODULEINIT    �M  \      SETUP    RO  H   a   SETUP%INIT    �O  H   a   SETUP%NEVAL "   �O  P   a   SETUP%DESCRIPTION    2P  P   a   SETUP%ATOMTYPE    �P  H   a   SETUP%NENV    �P  �   a   SETUP%ENVTYPES $   fQ  H   a   SETUP%NTYPES_GLOBAL    �Q  �   a   SETUP%GTYPE    BR  �   a   SETUP%LTYPE    �R  H   a   SETUP%RC_MIN    S  H   a   SETUP%RC_MAX    fS  P   a   SETUP%SFTYPE    �S  H   a   SETUP%NSF    �S  �   a   SETUP%SF    �T  H   a   SETUP%NSFPARAM    �T  �   a   SETUP%SFPARAM    �U  �   a   SETUP%SFENV     2V  �   a   SETUP%SFVAL_MIN     �V  �   a   SETUP%SFVAL_MAX     ZW  �   a   SETUP%SFVAL_AVG     �W  �   a   SETUP%SFVAL_COV    �X  �       SFVAL    Y  �       SFDERIV_I    �Y  �       SFDERIV_J    nZ  @       NSF_MAX    �Z  @       NNB_MAX 