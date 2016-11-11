#!/bin/bash

zonally_symmetric_tracer | lla2aeb q0a.vec
extract_field q0a.vec 0 | lla2aeb q0b.vec
correlate_fields q0a.vec q0b.vec

zonally_symmetric_tracer -E | lla2aeb q0a.vec
extract_field q0a.vec 0 | lla2aeb q0b.vec
correlate_fields q0a.vec q0b.vec

../libpetey/random 65160 | lla2aeb q0a.vec
extract_field q0a.vec 0 | lla2aeb q0b.vec
correlate_fields q0a.vec q0b.vec

rm q0a.vec q0b.vec

