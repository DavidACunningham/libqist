KPL/MK

   SPICE metakernel for basic trajectory use
   =====================================================================

   Original file name:                   mk_test_withorbit.tf
   Creation date:                        2024 Jun 25 11:41 CDT
   Created by:                           David Cunningham


		 Kernels to load are:

			Planetary Ephemeris SPK:      de440.bsp
			Moon rotating frame bpc:	  moon_pa_de440_200625.bpc
		    Moon rotating frame tf:		  moon_de440_220930.tf
			Test lunar orbit:			  test_orbit.bsp

			Leapseconds kernel (for
			time conversion):             naif0012.tls
		

	 \begindata
	 
	 PATH_VALUES = ('/home/david/wrk/nstgro/kernels'
				    '/home/david/wrk/nstgro/qist/kernels')

	 PATH_SYMBOLS = ('K'
					 'LK')

	 KERNELS_TO_LOAD = ( '$K/de440.bsp'
						 '$K/naif0012.tls'
						 '$K/moon_pa_de440_200625.bpc'
						 '$K/moon_de440_220930.tf'
						 '$LK/test_orbit.bsp'
						 )

	 \begintext

	 End of kernel
