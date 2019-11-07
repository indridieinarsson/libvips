/* local histogram equalisation
 *
 * Copyright: 1991, N. Dessipris 
 *
 * Author: N. Dessipris
 * Written on: 24/10/1991
 * Modified on : 
 * 25/1/96 JC
 *	- rewritten, adapting im_spcor()
 *	- correct result, 2x faster, partial, simpler, better arg checking
 * 8/7/04
 *	- expand input rather than output with new im_embed() mode
 *	- _raw() output is one pixel larger
 *	- sets Xoffset/Yoffset
 * 23/6/08
 * 	- check for window too small as well
 * 25/3/10
 * 	- gtkdoc
 * 	- small cleanups
 * 5/9/13
 * 	- redo as a class
 * 9/9/13
 * 	- any number of bands
 * 20/1/17
 * 	- add contrast limit
 * 	- sum to <= target, not < target, since cumulative hists include the
 * 	  current value
 * 	- scale result by 255, not 256, to avoid overflow
 * 	- off by 1 fix for odd window widths
 */

/*

    This file is part of VIPS.
    
    VIPS is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
    02110-1301  USA

 */

/*

    These files are distributed with VIPS - http://www.vips.ecs.soton.ac.uk

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/
#include <vips/intl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vips/vips.h>
#include <vips/internal.h>

typedef struct _VipsHistLocal {
	VipsOperation parent_instance;

	VipsImage *in;
	VipsImage *out;

	int width;
	int height;

	int max_slope;
    int n_bins;

} VipsHistLocal;

typedef VipsOperationClass VipsHistLocalClass;

G_DEFINE_TYPE( VipsHistLocal, vips_hist_local, VIPS_TYPE_OPERATION );

/* Our sequence value: the region this sequence is using, and local stats.
 */
typedef struct {
	VipsRegion *ir;		/* Input region */

	/* A 256-element hist for every band.
	 */
	unsigned int **hist;
} VipsHistLocalSequence;

static int
vips_hist_local_stop( void *vseq, void *a, void *b )
{
	VipsHistLocalSequence *seq = (VipsHistLocalSequence *) vseq;
	VipsImage *in = (VipsImage *) a;

	VIPS_UNREF( seq->ir );
	if( seq->hist &&
		in ) {
		int i; 

		for( i = 0; i < in->Bands; i++ )
			VIPS_FREE( seq->hist[i] );
		VIPS_FREE( seq->hist );
	}
	VIPS_FREE( seq );

	return( 0 );
}

static void *
vips_hist_local_start( VipsImage *out, void *a, void *b )
{
	VipsImage *in = (VipsImage *) a;
	const VipsHistLocal *local = (VipsHistLocal *) b;
	VipsHistLocalSequence *seq;

	int i;

	if( !(seq = VIPS_NEW( NULL, VipsHistLocalSequence )) )
		 return( NULL );
	seq->ir = NULL;
	seq->hist = NULL;

	if( !(seq->ir = vips_region_new( in )) || 
		!(seq->hist = VIPS_ARRAY( NULL, in->Bands, unsigned int * )) ) {
		vips_hist_local_stop( seq, NULL, NULL );
		return( NULL ); 
	}

	for( i = 0; i < in->Bands; i++ )
		if( !(seq->hist[i] = VIPS_ARRAY( NULL, local->n_bins, unsigned int )) ) {
		vips_hist_local_stop( seq, NULL, NULL );
		return( NULL ); 
	}

	return( seq );
}

/* Find histogram for the start of this line. 
 */
#define HISTOLOOP( ITYPE ) { \
    /*p1 = p; */ \
    printf("No Loop  done\n"); \
    int pelindex=0; \
    ITYPE *p = (ITYPE *) ppel; \
    ITYPE *q = (ITYPE *) qpel; \
    printf("Casting p1\n"); \
    VipsPel * ptmp = ppel; \
    ITYPE *p1; \
    p1 = (ITYPE *) ptmp; \
    printf("Done casting p1\n"); \
    for( j = 0; j < local->height; j++ ) { \
        for( i = 0, x = 0; x < local->width; x++ ) { \
            pelindex = p1[i]/npels_bin; \
            for( b = 0; b < bands; b++, i++ ) {\
                seq->hist[b][pelindex] += 1; \
            } \
        } \
        ptmp += lsk; \
        p1 = (ITYPE *) ptmp; \
    } \
    printf("Loop 1 done\n"); \
    /* Loop for output pels. */ \
    for( x = 0; x < r->width; x++ ){ \
        for( b = 0; b < bands; b++ ) { \
            /* Sum histogram up to current pel.  */ \
            unsigned int * hist = seq->hist[b]; \
            const int target = p[centre + b]; \
            /* This is the bin index of the target value in the histogram: */ \
            const int targetindex = target/npels_bin; \
            int sum; \
            sum = 0; \
            /* For CLAHE we need to limit the height of the */ \
            /* hist to limit the amount we boost the        */ \
            /* contrast by.                                 */ \
            if( max_slope > 0 ) { \
                int sum_over; \
                sum_over = 0; \
                /* Must be <= target, since a cum hist */ \
                /* always includes the current element.*/ \
                for( i = 0; i <= targetindex; i++ ) { \
                    if( hist[i] > max_slope ) { \
                        sum_over += hist[i] - max_slope; \
                        sum += max_slope; \
                    } \
                    else \
                    sum += hist[i]; \
                } \
                for( ; i < local->n_bins; i++ ) { \
                    if( hist[i] > max_slope ) \
                    sum_over += hist[i] - max_slope; \
                } \
                \
                /* The extra clipped off bit from the  */ \
                /* top of the hist is spread over all  */ \
                /* bins equally, then summed to target.*/ \
                /* ------------------------------------*/ \
                /* I'm not sure this is correct...     */ \
                sum += (target + 1) * sum_over / local->n_bins; \
            } \
            else { \
                sum = 0; \
                for( i = 0; i <= targetindex; i++ ) \
                    sum += hist[i]; \
            } \
            /* This can't overflow, even in        */ \
            /* contrast-limited mode.              */ \
            /* Scale by 255, not 256, or we'll get */ \
            /* overflow.                           */ \
            q[b] = maxpelval * sum / (local->width * local->height); \
            /* Adapt histogram --- remove the pels from */ \
            /* the left hand column, add in pels for a  */ \
            /* new right-hand column.                   */ \
            /* p1 = p + b; */ \
            ptmp = p + b; \
            p1 = (ITYPE *) ptmp; \
            for( j = 0; j < local->height; j++ ) { \
                hist[p1[0]] -= 1; \
                hist[p1[bands * local->width]] += 1; \
                ptmp += lsk; \
                p1 = (ITYPE *) ptmp; \
            } \
        } \
        p += bands; \
        q += bands; \
    } \
  }


static int
vips_hist_local_generate( VipsRegion *or, 
	void *vseq, void *a, void *b, gboolean *stop )
{
	VipsHistLocalSequence *seq = (VipsHistLocalSequence *) vseq;
	VipsImage *in = (VipsImage *) a;
	const VipsHistLocal *local = (VipsHistLocal *) b;
	VipsRect *r = &or->valid;
	const int bands = in->Bands; 
	const int max_slope = local->max_slope;
	const int n_bins = local->n_bins;

    printf("n_bins = %d\n", local->n_bins);
    printf("max_slope = %d\n", local->max_slope);

	VipsRect irect;
	int y;
	int lsk;
	int centre;		/* Offset to move to centre of window */
	/* What part of ir do we need?
	 */
    printf("debug 1 \n");
	irect.left = r->left;
	irect.top = r->top;
	irect.width = r->width + local->width; 
	irect.height = r->height + local->height; 
    printf("debug 2 %d, %d\n", irect.width, irect.height );
	if( vips_region_prepare( seq->ir, &irect ) ) {
        printf("prepare failed : \n");
		return( -1 );
    }

    printf("debug 3 \n");
	lsk = VIPS_REGION_LSKIP( seq->ir );
    printf("debug 4 \n");
	centre = lsk * (local->height / 2) + bands * (local->width / 2);
    printf("debug 5 \n");

	for( y = 0; y < r->height; y++ ) {
		/* Get input and output pointers for this line.
		 */
        printf("Get ppel\n");
		VipsPel * ppel = 
			VIPS_REGION_ADDR( seq->ir, r->left, r->top + y );
        printf("Get qpel\n");
		VipsPel * qpel = 
			VIPS_REGION_ADDR( or, r->left, r->top + y );
		// VipsPel * restrict p1;

		int x, i, j, b;

		/* Find histogram for the start of this line. 
		 */
        printf("zero bins ppel\n");
		for( b = 0; b < bands; b++ )
			memset( seq->hist[b], 0, local->n_bins * sizeof( unsigned int ) );

        int maxpelval = USHRT_MAX;
        //int maxpelval = UCHAR_MAX;
        printf("maxpelval = %d\n", maxpelval);
        printf("n_bins = %d\n", local->n_bins);
        int npels_bin = maxpelval/local->n_bins;
        printf("npels_bin = %d\n", npels_bin);
        printf("Call HISTOLOOP\n");
        HISTOLOOP(unsigned short);
        //HISTOLOOP(unsigned char);
	}
	return( 0 );
}

static int
vips_hist_local_build( VipsObject *object )
{
	VipsObjectClass *class = VIPS_OBJECT_GET_CLASS( object );
	VipsHistLocal *local = (VipsHistLocal *) object;
	VipsImage **t = (VipsImage **) vips_object_local_array( object, 3 );

	VipsImage *in;

	if( VIPS_OBJECT_CLASS( vips_hist_local_parent_class )->build( object ) )
		return( -1 );

	in = local->in; 

	if( vips_image_decode( in, &t[0] ) )
		return( -1 );
	in = t[0]; 

    switch( vips_image_get_format( in ) ) {
        case VIPS_FORMAT_UCHAR:
            printf("uchar : %d\n", 256);
            break;
        case VIPS_FORMAT_USHORT:
            printf("ushort : %d\n", 65536);
            break;
    }
  //local->npels_bin = local->maxpelval/local->n_bins;
  //printf("Setting npels_bin %d, %d, %d\n",local->maxpelval, local->n_bins, local->npels_bin);

  //if( vips_check_format( class->nickname, in, VIPS_FORMAT_UCHAR ) )
  //if( vips_check_format( class->nickname, in, VIPS_FORMAT_USHORT ) )
  //  return( -1 );

	if( local->width > in->Xsize || 
		local->height > in->Ysize ) {
		vips_error( class->nickname, "%s", _( "window too large" ) );
		return( -1 );
	}

	/* Expand the input. 
	 */
  printf("width %d x height %d \n", local->width , local->height );
  printf("Other arguments : %d, %d, %d, %d\n", in->Xsize, local->width , in->Ysize , local->height);
	if( vips_embed( in, &t[1], 
		local->width / 2, local->height / 2, 
		in->Xsize + local->width - 1, in->Ysize + local->height - 1,
                  //"extend", VIPS_EXTEND_WHITE,
    "extend", VIPS_EXTEND_MIRROR,
    NULL ) )
    return( -1 );
	in = t[1];

	g_object_set( object, "out", vips_image_new(), NULL ); 

	/* Set demand hints. FATSTRIP is good for us, as THINSTRIP will cause
	 * too many recalculations on overlaps.
	 */
  //  VIPS_DEMAND_STYLE_FATSTRIP, in, NULL ) )
	if( vips_image_pipelinev( local->out, 
		VIPS_DEMAND_STYLE_FATSTRIP, in, NULL ) )
		return( -1 );
	local->out->Xsize -= local->width - 1;
	local->out->Ysize -= local->height - 1;

	if( vips_image_generate( local->out, 
		vips_hist_local_start, 
		vips_hist_local_generate, 
		vips_hist_local_stop, 
		in, local ) )
		return( -1 );

	local->out->Xoffset = 0;
	local->out->Yoffset = 0;

	vips_reorder_margin_hint( local->out, local->width * local->height ); 

	return( 0 );
}

static void
vips_hist_local_class_init( VipsHistLocalClass *class )
{
	GObjectClass *gobject_class = G_OBJECT_CLASS( class );
	VipsObjectClass *object_class = (VipsObjectClass *) class;

	gobject_class->set_property = vips_object_set_property;
	gobject_class->get_property = vips_object_get_property;

	object_class->nickname = "hist_local";
	object_class->description = _( "local histogram equalisation" );
	object_class->build = vips_hist_local_build;

	VIPS_ARG_IMAGE( class, "in", 1, 
		_( "Input" ), 
		_( "Input image" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsHistLocal, in ) );

	VIPS_ARG_IMAGE( class, "out", 2, 
		_( "Output" ), 
		_( "Output image" ),
		VIPS_ARGUMENT_REQUIRED_OUTPUT, 
		G_STRUCT_OFFSET( VipsHistLocal, out ) );

	VIPS_ARG_INT( class, "width", 4, 
		_( "Width" ), 
		_( "Window width in pixels" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsHistLocal, width ),
		1, VIPS_MAX_COORD, 1 );

	VIPS_ARG_INT( class, "height", 5, 
		_( "Height" ), 
		_( "Window height in pixels" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsHistLocal, height ),
		1, VIPS_MAX_COORD, 1 );

	VIPS_ARG_INT( class, "max_slope", 6, 
		_( "Max slope" ), 
		_( "Maximum slope (CLAHE)" ),
		VIPS_ARGUMENT_OPTIONAL_INPUT,
		G_STRUCT_OFFSET( VipsHistLocal, max_slope ),
		0, 100, 0 );

	VIPS_ARG_INT( class, "n_bins", 7, 
                _( "Number of bins" ), 
                _( "Number of histogram bins" ),
                VIPS_ARGUMENT_OPTIONAL_INPUT,
                G_STRUCT_OFFSET( VipsHistLocal, n_bins ),
                0, 65536, 256 );
}

static void
vips_hist_local_init( VipsHistLocal *local )
{
}

/**
 * vips_hist_local: (method)
 * @in: input image
 * @out: (out): output image
 * @width: width of region
 * @height: height of region
 * @...: %NULL-terminated list of optional named arguments
 *
 * Optional arguments:
 *
 * * @max_slope: maximum brightening
 * * @n_bins: Number of histogram bins
 *
 * Performs local histogram equalisation on @in using a
 * window of size @width by @height centered on the input pixel. 
 *
 * The output image is the same size as the input image. The edge pixels are
 * created by mirroring the input image outwards.
 *
 * If @max_slope is greater than 0, it sets the maximum value for the slope of
 * the cumulative histogram, that is, the maximum brightening that is
 * performed. A value of 3 is often used. Local histogram equalization with
 * contrast limiting is usually called CLAHE.
 *
 * See also: vips_hist_equal().
 *
 * Returns: 0 on success, -1 on error
 */
int 
vips_hist_local( VipsImage *in, VipsImage **out, int width, int height, ... )
{
	va_list ap;
	int result;

	va_start( ap, height );
	result = vips_call_split( "hist_local", ap, in, out, width, height );
	va_end( ap );

	return( result );
}
