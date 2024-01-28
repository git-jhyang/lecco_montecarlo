#!/bin/sh

diff sample.output check.output
if [ $? -ne 0 ]; then
  echo "f_jump_coeff produced incorrect results."
  echo "Check Error."
  exit -1
else
  echo "f_jump_coeff produced correct results."
  echo "Check OK."
  exit 0
fi

