/****************************/
/* THIS IS OPEN SOURCE CODE */
/****************************/

/** 
 * @author  Vince Weaver
 *
 * test case for micpower component 
 * Based on coretemp test code by Vince Weaver
 * 
 *
 * @brief
 *   Tests basic component functionality
 */

#include <stdio.h>
#include <stdlib.h>
#include "papi_test.h"

#define NUM_EVENTS 1
long long g_time_e;
long long g_time_s;

int RCRMICPowerCheck()
{

        int retval,cid,numcmp;
	int EventSet = PAPI_NULL;
	long long values[NUM_EVENTS];
	int code;
	char event_name[PAPI_MAX_STR_LEN];
	int total_events=0;
	int r;
	const PAPI_component_info_t *cmpinfo = NULL;

    ///* Set TESTS_QUIET variable */
    //tests_quiet( argc, argv );      
    g_time_s = PAPI_get_real_usec();

	/* PAPI Initialization */
	retval = PAPI_library_init( PAPI_VER_CURRENT );
	if ( retval != PAPI_VER_CURRENT ) {
	   test_fail(__FILE__, __LINE__,"PAPI_library_init failed\n",retval);
	}

        numcmp = PAPI_num_components();

	for(cid=0; cid<numcmp; cid++) {

	   if (!TESTS_QUIET) {
	      if ( (cmpinfo = PAPI_get_component_info(cid)) == NULL) {
	         test_fail(__FILE__, __LINE__,"PAPI_get_component_info failed\n", 0);
	      }
	      //printf("\tComponent %d - %s\n", cid, cmpinfo->name);
	   }

       if ( 0 != strncmp(cmpinfo->name,"micpower",8)) {
            continue;
       }

	   code = PAPI_NATIVE_MASK;
      while (1) {
           r = PAPI_enum_cmp_event( &code, PAPI_ENUM_FIRST, cid );

	   while ( r == PAPI_OK ) {
	      retval = PAPI_event_code_to_name( code, event_name );
	      if ( retval != PAPI_OK ) {
		 printf("Error translating %#x\n",code);
	         test_fail( __FILE__, __LINE__, 
                            "PAPI_event_code_to_name", retval );
	      }
          if (0 != strncmp(event_name, "micpower:::tot0", 15) ) {
                r = PAPI_enum_cmp_event( &code, PAPI_ENUM_EVENTS, cid );
                continue;
          }

	      //if (!strncmp(event_name,"micpower:::tot0", 15)) {
	         if (!TESTS_QUIET) printf("%#x %s ",code,event_name);
	     
	         EventSet = PAPI_NULL;

	         retval = PAPI_create_eventset( &EventSet );
	         if (retval != PAPI_OK) {
	            test_fail(__FILE__, __LINE__, 
                              "PAPI_create_eventset()",retval);
		 }

	         retval = PAPI_add_event( EventSet, code );
	         if (retval != PAPI_OK) {
	            test_fail(__FILE__, __LINE__, 
                                 "PAPI_add_event()",retval);
		 }

	         retval = PAPI_start( EventSet);
	         if (retval != PAPI_OK) {
	            test_fail(__FILE__, __LINE__, "PAPI_start()",retval);
		 }

	         retval = PAPI_stop( EventSet, values);
	         if (retval != PAPI_OK) {
	            test_fail(__FILE__, __LINE__, "PAPI_start()",retval);
		 }

	         if (!TESTS_QUIET) printf(" value: %lld\n",values[0]);
             g_time_e = PAPI_get_real_usec();
             // store to shared memory
             update_blackboard(event_name, values[0], g_time_e-g_time_s);
             g_time_s = g_time_e;

	         retval = PAPI_cleanup_eventset( EventSet );
	         if (retval != PAPI_OK) {
	            test_fail(__FILE__, __LINE__, 
                              "PAPI_cleanup_eventset()",retval);
		 }

	         retval = PAPI_destroy_eventset( &EventSet );
	         if (retval != PAPI_OK) {
	            test_fail(__FILE__, __LINE__, 
                              "PAPI_destroy_eventset()",retval);
		 }

	         total_events++;
	      //}
	      r = PAPI_enum_cmp_event( &code, PAPI_ENUM_EVENTS, cid );
	   }
       usleep(50000);
      }
     }

	if (total_events==0) {

	   test_skip(__FILE__,__LINE__,"No coretemp events found",0);
	}

//	test_pass( __FILE__, NULL, 0 );
		
	return 0;
}

