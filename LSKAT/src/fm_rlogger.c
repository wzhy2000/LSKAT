/* logger.c  -	log functions
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

#include <R.h>
#include <Rinternals.h>

#include <Rembedded.h>
#include <Rdefines.h>

#include "fm_rlogger.h"


//typedef unsigned char BYTE;

static char szPid[256];
static char* pgm_name  = NULL;
static char* szLogName = NULL;
static int  nErrorCount= 0;
static FILE* pFLog     = NULL;
static bool bDebug     = false;

/****************
 * start the log file, if szLogFile is NULL, stdout or stderr is used to output log.
 * the Fd where logoutputs should go.
 */

int start_log(bool bCmdDebug)
{
    bDebug = bCmdDebug;
    if( bDebug )
       	pFLog = stdout;
	else
   		pFLog = stderr;

    return(0);
}

int stop_log()
{
    return(0);
}

void log_set_pid( int pid )
{
    if( pid )
        sprintf(szPid,"[%u]", (unsigned)pid );
    else
        *szPid = 0;
}

FILE* get_log_stream()
{
    if( !pFLog )
        pFLog = stderr;
    return pFLog;
}

const char* get_log_name(void)
{
    return szLogName? szLogName : "";
}

int get_log_errorcount( int clear)
{
    int n = nErrorCount;
    if( clear )
        nErrorCount = 0;

    return n;
}

void log_print_prefix(const char *text)
{
    if( pgm_name )
        fprintf(pFLog, "%s%s: %s", pgm_name, szPid, text );
    else
        fprintf(pFLog, "?%s: %s", szPid, text );
}

static void log_print_prefix_f(const char *text, const char *fname)
{
    if( pgm_name )
        fprintf(pFLog, "%s%s:%s: %s", pgm_name, szPid, fname, text );
    else
        fprintf(pFLog, "?%s:%s: %s", szPid, fname, text );
}

void _log_debug( const char* szSrc, int nSrcLine, const char* fmt, ... )
{
    va_list arg_ptr ;
    if (bDebug)
    {
        log_print_prefix_f("DBG: ", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
        fprintf(pFLog, "\n");
        fflush(pFLog);
    }
}

void _log_info( const char* szSrc, int nSrcLine, const char* fmt, ... )
{
    va_list arg_ptr ;
    va_start( arg_ptr, fmt ) ;
    vfprintf( stderr, fmt,arg_ptr) ;
    fprintf( stderr, "\n");
    va_end(arg_ptr);

/*
    if (bDebug)
    {
        va_list arg_ptr ;
        log_print_prefix_f("", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
    }

    fprintf(pFLog, "\n");
    fflush(pFLog);
*/
}


void _log_prompt( const char* szSrc, int nSrcLine, const char* fmt, ... )
{
    //output the contens to console
    va_list arg_ptr ;
    va_start( arg_ptr, fmt ) ;
    vfprintf( stderr, fmt,arg_ptr) ;
    fprintf( stderr, "\n");
    va_end(arg_ptr);

/*    if (bDebug)
    {
        va_list arg_ptr ;
        log_print_prefix_f("", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
    }
    else
    {
        va_list arg_ptr ;
        log_print_prefix("");

        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
    }

    fprintf(pFLog, "\n");
    fflush(pFLog);
*/

}

void _log_error( const char* szSrc, int nSrcLine, const char*  fmt, ... )
{
    va_list arg_ptr ;
    va_start( arg_ptr, fmt ) ;
    vfprintf( stderr, fmt,arg_ptr) ;
    fprintf( stderr, "\n");
    va_end(arg_ptr);

/*    if (bDebug)
    {
        va_list arg_ptr ;
        log_print_prefix_f("", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
        nErrorCount++;
    }
    else
    {
        va_list arg_ptr ;
        log_print_prefix("");
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
        nErrorCount++;
    }

    fprintf(pFLog, "\n");
    fflush(pFLog);
*/
}

void _log_fatal( const char* szSrc, int nSrcLine, const char*  fmt, ... )
{
    va_list arg_ptr ;
    va_start( arg_ptr, fmt ) ;
    vfprintf( stderr, fmt,arg_ptr) ;
    fprintf( stderr, "\n");
    va_end(arg_ptr);

/*    if (bDebug)
    {
        va_list arg_ptr ;
        log_print_prefix_f("fatal: ", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
        //secmem_dump_stats();
    }
    else
    {
        va_list arg_ptr ;
        log_print_prefix("fatal: ");
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
        //secmem_dump_stats();
    }

    fprintf(pFLog, "\n");
    fflush(pFLog);
*/

    stop_log();
    error("An Exception was caugth, please see log file.");
}

void _log_hexdump( const char *text, const char *buf, size_t len )
{
    int i;

    log_print_prefix(text);
    for(i=0; i < len; i++ )
    {
        char c = buf[i];
        fprintf(pFLog, " %02X", (unsigned int)c );
    }

    fputc('\n', pFLog);
}
