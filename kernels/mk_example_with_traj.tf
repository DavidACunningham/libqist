KPL/TF

   SPICE metakernel for QIST example
   =====================================================================

   Original file name:                   mk_example_with_traj.tf
   Creation date:                        2024 Nov 22 12:20 CST
   Created by:                           David Cunningham


		 Kernels to load are:

			Planetary Ephemeris SPK:      de440.bsp
			Mars System Kernel:			  mar097.bsp
			Deimos orbiter trajectory:    example_binary_kernel.bsp

			Leapseconds kernel (for
			time conversion):             naif0012.tls
			Deimos Rotating Frame Data:   pck00011_n0066.tpc
			Mars-Deimos Rotating Frame:   mdrot.tf

	 \begindata
	 
	 PATH_VALUES = ('C:\Users\tester\Downloads\libqist\libqist\kernels')

	 PATH_SYMBOLS = ('K')

	 KERNELS_TO_LOAD = ( '$K/de440.bsp'
                       '$K/mar097.bsp'
                       '$K/example_binary_kernel.bsp'
                       '$K/naif0012.tls.pc'
                       '$K/mdrot.tf'
                       )	 

	 \begintext

	 End of kernel
