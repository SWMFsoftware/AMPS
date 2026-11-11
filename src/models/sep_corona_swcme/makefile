CXX ?= g++
AR ?= ar
PYTHON ?= python3
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic -Werror
# ../../general/input supplies the header-only generic AMPS runtime-input
# reader (amps_input.h; grammar in src/INPUT_FILE.md).  It has no PIC, MPI or
# application dependency, so the dependency-light layering is preserved.
CPPFLAGS := -Iinclude -Ishock_front -I../sep_coronal_cme/include \
	-I../sep_coronal_cme/src -I../sep_common -I../../general/input
BUILD := build
SOURCES := $(wildcard src/*.cpp)
SHOCK_FRONT_SOURCES := $(wildcard shock_front/*.cpp)
OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SOURCES)) \
	$(patsubst shock_front/%.cpp,$(BUILD)/shock_front_%.o,$(SHOCK_FRONT_SOURCES))
CORONAL_ARCHIVE := ../sep_coronal_cme/build/libsep_coronal_cme.a

.PHONY: all lib test clean coronal-lib check-architecture shock-front-test
.DEFAULT_GOAL := all

all: lib

lib: $(BUILD)/libsep_corona_swcme.a

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/shock_front_%.o: shock_front/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/libsep_corona_swcme.a: $(OBJECTS)
	$(AR) rcs $@ $^

coronal-lib:
	$(MAKE) -C ../sep_coronal_cme -j16 lib

$(BUILD)/test_bg3d1: test/test_bg3d1.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d2: test/test_bg3d2.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d3: test/test_bg3d3.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4: test/test_bg3d4.cpp test/bg3d4_reference_fixtures.h \
        test/contact_derivative_reference.h \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_piston_contact: test/test_bg3d4_piston_contact.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_planar_piston: test/test_bg3d4_planar_piston.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_planar_formation: test/test_bg3d4_planar_formation.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_planar_mhd: test/test_bg3d4_planar_mhd.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_spherical_piston: test/test_bg3d4_spherical_piston.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_ambient_balance: test/test_bg3d4_ambient_balance.cpp \
		test/shell_quadrature_reference.h \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_append: test/test_bg3d4_append.cpp \
		test/shell_quadrature_reference.h \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_piston_sheath: test/test_bg3d4_piston_sheath.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4_piston_refinement: test/test_bg3d4_piston_refinement.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

# This test shares only the independently differentiated-position grader with
# CMBGU05. Its manufactured trajectories and imposed slip are separate from
# the production material map, so agreement cannot be an API self-comparison.
$(BUILD)/test_bg3d4_contact_derivative: test/test_bg3d4_contact_derivative.cpp \
        test/contact_derivative_reference.h coronal-lib | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_spherical_shell: test/test_spherical_shell.cpp \
        include/sep_corona_swcme/spherical_shell.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(BUILD)/test_cross_ray_transport: test/test_cross_ray_transport.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_spherical_sector: test/test_spherical_sector.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_cross_ray_mesh: test/test_cross_ray_mesh.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_cross_ray_ct: test/test_cross_ray_ct.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_cross_ray_reconstruction: test/test_cross_ray_reconstruction.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_cross_ray_fluid: test/test_cross_ray_fluid.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_cross_ray_contact: test/test_cross_ray_contact.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_material_transport: test/test_material_transport.cpp \
        $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_shock_front: test/shock_front/test_reduced.cpp \
		$(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

shock-front-test: $(BUILD)/test_shock_front
	$(BUILD)/test_shock_front

check-architecture: $(BUILD)/libsep_corona_swcme.a
	$(PYTHON) test/check_architecture.py

test: $(BUILD)/test_bg3d1 $(BUILD)/test_bg3d2 $(BUILD)/test_bg3d3 \
		$(BUILD)/test_bg3d4 $(BUILD)/test_bg3d4_piston_contact \
		$(BUILD)/test_bg3d4_planar_piston \
		$(BUILD)/test_bg3d4_planar_formation \
		$(BUILD)/test_bg3d4_planar_mhd \
		$(BUILD)/test_bg3d4_spherical_piston \
		$(BUILD)/test_bg3d4_ambient_balance \
		$(BUILD)/test_bg3d4_append \
		$(BUILD)/test_bg3d4_piston_sheath \
		$(BUILD)/test_bg3d4_piston_refinement \
		$(BUILD)/test_bg3d4_contact_derivative \
		$(BUILD)/test_spherical_shell $(BUILD)/test_cross_ray_transport \
		$(BUILD)/test_spherical_sector $(BUILD)/test_cross_ray_mesh \
		$(BUILD)/test_cross_ray_ct $(BUILD)/test_cross_ray_reconstruction \
		$(BUILD)/test_cross_ray_fluid $(BUILD)/test_cross_ray_contact \
		$(BUILD)/test_material_transport check-architecture
	$(BUILD)/test_bg3d1
	$(BUILD)/test_bg3d2
	$(BUILD)/test_bg3d3
	$(BUILD)/test_bg3d4
	$(BUILD)/test_bg3d4_piston_contact
	$(BUILD)/test_bg3d4_planar_piston
	$(BUILD)/test_bg3d4_planar_formation
	$(BUILD)/test_bg3d4_planar_mhd
	$(BUILD)/test_bg3d4_spherical_piston
	$(BUILD)/test_bg3d4_ambient_balance
	$(BUILD)/test_bg3d4_append
	$(BUILD)/test_bg3d4_piston_sheath
	$(BUILD)/test_bg3d4_piston_refinement
	$(BUILD)/test_bg3d4_contact_derivative
	$(BUILD)/test_spherical_shell
	$(BUILD)/test_cross_ray_transport
	$(BUILD)/test_spherical_sector
	$(BUILD)/test_cross_ray_mesh
	$(BUILD)/test_cross_ray_ct
	$(BUILD)/test_cross_ray_reconstruction
	$(BUILD)/test_cross_ray_fluid
	$(BUILD)/test_cross_ray_contact
	$(BUILD)/test_material_transport

clean:
	rm -rf $(BUILD)

-include $(OBJECTS:.o=.d)
