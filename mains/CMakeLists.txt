ecbuild_add_executable( TARGET  bump.x
                        SOURCES bump.F90
                                type_model.F90
                                type_timer.F90
                        LIBS    bump
                       )
