KPL/TF

   SPICE metakernel for QIST example
   =====================================================================

   Original file name:                   mk_example.tf
   Creation date:                        2024 Nov 22 12:20 CST
   Created by:                           David Cunningham


		 Kernels to load are:

			Planetary Ephemeris SPK:      de440.bsp
			Mars System Kernel:			  mar097.bsp
			Planetary Constants Kernel:   pck00011.tpc

			Leapseconds kernel (for
			time conversion):             naif0012.tls
			Mars-Deimos Rotating Frame:   mdrot.tf

	 \begindata
	 
	 PATH_VALUES = ('/home/david/wrk/nstgro/libqist/kernels')

	 PATH_SYMBOLS = ('K')

	 KERNELS_TO_LOAD = ( '$K/de440.bsp'
						 '$K/mar097.bsp'
						 '$K/naif0012.tls'
						 '$K/pck00011.tpc'
						 '$K/mdrot.tf'
					   )	 

	 \begintext

	 End of kernel
