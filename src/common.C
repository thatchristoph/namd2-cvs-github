/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   global functions as declared in common.h
*/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>

#include "common.h"
#include "InfoStream.h"

#include "charm++.h"

#ifdef WIN32
#include <io.h>
#define access(PATH,MODE) _access(PATH,00)
#define NOCOMPRESSED
#endif

#ifdef USESTDIOSTUBS

// On the ASCI-Red, these functions have prototypes that don't match
// the standard, so I'll override them

int write(int fd, const void* buf, int count) 
{ 
  return write(fd,(void*)buf,(size_t)count);
}

int stat(const char* fn, struct stat* buf)
{
  return stat((char*)fn, buf);
}

#endif


// make a duplicate of a string
char *NAMD_stringdup(const char *s) {
  char *rs;

  if(!s)
    return NULL;

  rs = new char[strlen(s) + 1];
  strcpy(rs,s);

  return rs;
}


// signal all nodes, it's time to quit
void NAMD_quit(const char *err_msg)

{
   char *new_err_msg = new char[strlen(err_msg) + 20];
   sprintf(new_err_msg,"EXITING: %s\n",err_msg);
   CkPrintf(new_err_msg);
   fflush(stdout);
   CmiAbort(new_err_msg);
   delete [] new_err_msg;
}

 
// signal all nodes, it's time to quit
void NAMD_die(const char *err_msg)

{
   char *new_err_msg = new char[strlen(err_msg) + 20];
   sprintf(new_err_msg,"FATAL ERROR: %s\n",err_msg);
   CkPrintf(new_err_msg);
   fflush(stdout);
   CmiAbort(new_err_msg);
   delete [] new_err_msg;
}


// signal all nodes, it's time to quit
void NAMD_err(const char *err_msg)

{
   char *sys_err_msg = strerror(errno);
   if ( ! sys_err_msg ) sys_err_msg = "(unknown error)";
   char *new_err_msg = new char[strlen(err_msg) + 20 + strlen(sys_err_msg)];
   sprintf(new_err_msg,"FATAL ERROR: %s: %s\n",err_msg, sys_err_msg);
   CkPrintf(new_err_msg);
   fflush(stdout);
   CmiAbort(new_err_msg);
   delete [] new_err_msg;
}


// signal all nodes, it's time to quit and it's our own damn fault
void NAMD_bug(const char *err_msg)

{
   const char *bug_msg = 
     "FATAL ERROR: See http://www.ks.uiuc.edu/Research/namd/bugreport.html";
   char *new_err_msg = new char[strlen(err_msg) + 20 + strlen(bug_msg)];
   sprintf(new_err_msg,"FATAL ERROR: %s\n%s\n",err_msg,bug_msg);
   CkPrintf(new_err_msg);
   fflush(stdout);
   CmiAbort(new_err_msg);
   delete [] new_err_msg;
}

// move filename to filename.BAK
void NAMD_backup_file(const char *filename, const char *extension)
{
  if (access(filename, F_OK) == 0) {
    if ( ! extension ) extension = ".BAK";
    char *backup = new char[strlen(filename)+strlen(extension)+1];
    strcpy(backup, filename);
    strcat(backup, extension);
#if defined(WIN32) && !defined(__CYGWIN__)
    if ( remove(backup) ) if ( errno != ENOENT ) {
      char *sys_err_msg = strerror(errno);
      if ( ! sys_err_msg ) sys_err_msg = "(unknown error)";
      iout << iERROR << "Error on removing file "
	<< backup << ": " << sys_err_msg << "\n" << endi;
      fflush(stdout);
    }
#endif
    if ( rename(filename,backup) )
    {
      char *sys_err_msg = strerror(errno);
      if ( ! sys_err_msg ) sys_err_msg = "(unknown error)";
      iout << iERROR << "Error on renaming file " << filename
	<< " to " << backup << ": " << sys_err_msg << "\n" << endi;
      fflush(stdout);
      // char errmsg[256];
      // sprintf(errmsg, "Error on renaming file %s to %s",filename,backup);
      // NAMD_err(errmsg);
    }
    delete [] backup;
  }
}

// same as write, only does error checking internally
void NAMD_write(int fd, const char *buf, size_t count) {
  while ( count ) {
#if defined(WIN32) && !defined(__CYGWIN__)
    long retval = _write(fd,buf,count);
#else
    ssize_t retval = write(fd,buf,count);
#endif
    if ( retval < 0 && errno == EINTR ) retval = 0;
    if ( retval < 0 ) NAMD_die(strerror(errno));
    if ( retval > count ) NAMD_bug("extra bytes written in NAMD_write()");
    buf += retval;
    count -= retval;
  }
}


/***************************************************************************
 Fopen(char *Filename, char *mode):  similar to fopen(filename,mode) except
 it checks for compressed file names too.
 For example:  Fopen("config");
   This will first look for the filename "config" (and return "r" file handle
   if it is found).
   Then it will look for "config.Z" (and run "zcat config.Z", returning
   a file handle to the uncompressed data if found).
   Then it will look for "config.gz" (and run "gzip -d -c config.gz", returning
   a file handle to the uncompressed data if found).
 ***************************************************************************/
FILE *Fopen	(const char *filename, const char *mode)
{
  struct stat buf;
  // check if basic filename exists (and not a directory)

#if defined(NOCOMPRESSED)
  if (!stat(filename,&buf))
    {
      FILE *rval;
      while ( ! (rval = fopen(filename,mode)) ) {
        if ( errno != EINTR ) break;
      }
      return(rval);
    }
#else
  if (!stat(filename,&buf))
    {
      if (!S_ISDIR(buf.st_mode)) {
        FILE *rval;
        while ( ! (rval = fopen(filename,mode)) ) {
          if ( errno != EINTR ) break;
        }
        return(rval);
      }
    }
  // check for a compressed file
  char *realfilename;
  char *command;
  FILE *fout;
  command = (char *)malloc(strlen(filename)+25);
  // check for .Z (unix compress)
  sprintf(command,"zcat %s.Z",filename);
  realfilename = command+5;
  iout << "Command = " << command << "\n" << endi;
  iout << "Filename.Z = " << realfilename << "\n" << endi;
  if (!stat(realfilename,&buf))
	{
	if (!S_ISDIR(buf.st_mode))
		{
		fout = popen(command,mode);
		// on HP-UX, the first character(s) out of pipe may be
		// garbage!  (Argh!)
		int C;
		do
		  {
		  C = fgetc(fout);
		  // iout << "C is " << C << "\n" << endi;
		  if (isalnum(C) || isspace(C))
			{
			ungetc(C,fout);
			C = -1;	// outta loop
			}
		  } while(C != -1);
		free(command);
		return(fout);
		}
	}
  // check for .gz (gzip)
  sprintf(command,"gzip -d -c %s.gz",filename);
  realfilename = command+11;
  iout << "Command = " << command << "\n" << endi;
  iout << "Filename.gz = " << realfilename << "\n" << endi;
  if (!stat(realfilename,&buf))
	{
	if (!S_ISDIR(buf.st_mode))
		{
		fout = popen(command,mode);
		// on HP-UX, the first character(s) out of pipe may be
		// garbage!  (Argh!)
		int C;
		do
		  {
		  C = fgetc(fout);
		  // iout << "C is " << C << "\n" << endi;
		  if (isalnum(C) || isspace(C))
			{
			ungetc(C,fout);
			C = -1;	// outta loop
			}
		  } while(C != -1);
		free(command);
		return(fout);
		}
	}
  free(command);
#endif /* !defined(NOCOMPRESSED) */

  return(NULL);
} /* Fopen() */

/***************************************************************************
 Fclose(FILE *fout):  similar to fclose(fout) except it first checks if the
 file handle fout is a named pipe.
 ***************************************************************************/
int	Fclose	(FILE *fout)
{
  int rc = -1;
#if !defined(NOCOMPRESSED)
  rc = pclose(fout);
#endif
  if (rc == -1)	// stream not associated with a popen()
    {
    rc = fclose(fout);
    }
  return rc;
} /* Fclose() */


