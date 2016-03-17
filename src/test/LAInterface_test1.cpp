#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"
#include "LAInterface.h"

#include "gtest/gtest.h"

#include <iostream>

TEST(SimpleTest, Plus) {
	ASSERT_EQ(1+1, 2);
}