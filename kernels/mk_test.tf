KPL/FK

   SPICE metakernel for basic trajectory use
   =====================================================================

   Original file name:                   mk.tf
   Creation date:                        2023 Apr 13 17:06 CDT
   Created by:                           David Cunningham


		 Kernels to load are:

			Planetary Ephemeris SPK:      de440.bsp
			Planetary Constants Kernel:   pck00011.tpc
			Mass Parameters Kernel:       de-403-masses.tpc
			Moon rotating frame bpc:	  moon_pa_de440_200625.bpc
		    Moon rotating frame tf:		  moon_de440_220930.tf

			Leapseconds kernel (for
			time conversion):          naif0012.tls
		
			Earth-Moon Rotating Frame: emrot.tf
			Spacecraft Kernel:		perturbed_ker.bsp

	 \begindata
	 
	 PATH_VALUES = ('/home/david/wrk/nstgro/kernels'
				    '/home/david/wrk/nstgro/qist/kernels')

	 PATH_SYMBOLS = ('K'
					 'LK')

	 KERNELS_TO_LOAD = ( '$K/de440.bsp'
						 '$K/naif0012.tls'
						 '$K/pck00011.tpc'
						 '$LK/twobod_ker.bsp'
						 '$K/moon_pa_de440_200625.bpc'
						 '$K/moon_de440_220930.tf'
					   )	 

	 \begintext

	 End of kernel
