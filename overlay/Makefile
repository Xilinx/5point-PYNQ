DEVICE := xilinx_u200_xdma_201830_2
TARGET := hw
FREQUENCY := 170
VPP := v++

CP = cp -rf
RM = rm -f
RMDIR = rm -rf

# device2xsa - create a filesystem friendly name from device name
# $(1) - full name of device
device2xsa = $(strip $(patsubst %.xpfm, % , $(shell basename $(DEVICE))))
XSA := $(call device2xsa, $(DEVICE))

XO_DIR := ./_xo.$(TARGET).$(XSA)
XCLBIN_DIR := ./_xclbin.$(TARGET).$(XSA)
SRC_DIR := ./src
INSTALL_PATH := ../fivepoint_pynq/notebooks

KERNEL_XO := $(XO_DIR)/5point_nister.xo
KERNEL_NAME := fivept_nister

INSTALL_TARGETS := $(INSTALL_PATH)/server/5point_nister.xclbin

# Kernel compiler global settings
CLFLAGS += -O3 -t $(TARGET) --platform $(DEVICE)
ifneq ($(TARGET), hw)
	CLFLAGS += -g
endif

# LDCLFLAGS += --hls.max_memory_ports $(KERNEL_NAME)  # valid only for OpenCL kernels
# LDCLFLAGS += --hls.memory_port_data_width $(KERNEL_NAME):64  # valid only for OpenCL kernels
# LDCLFLAGS += --connectivity.sp $(KERNEL_NAME)_1:m_axi_gmem0:bank0  # depend on 'hls.max_memory_ports'
# LDCLFLAGS += --connectivity.sp $(KERNEL_NAME)_1:m_axi_gmem1:bank1  # depend on 'hls.max_memory_ports'
# LDCLFLAGS += --connectivity.sp $(KERNEL_NAME)_1:m_axi_gmem2:bank2  # depend on 'hls.max_memory_ports'

all: build install

build: check-vitis check-xrt 5point_nister.$(TARGET).$(XSA).xclbin

# install targets assume both xclbin and dst folder exist
install: $(INSTALL_TARGETS)

check-vitis:
ifndef XILINX_VITIS
	$(error XILINX_VITIS is not set. Please make sure you have sourced the Vitis settings64.{csh,sh})
endif

check-xrt:
ifndef XILINX_XRT
	$(error XILINX_XRT variable is not set. Please make sure you have sourced the XRT setup.{csh,sh})
endif

$(INSTALL_PATH)/server/5point_nister.xclbin:
ifneq (,$(wildcard 5point_nister.$(TARGET).$(XSA).xclbin))
ifneq (,$(wildcard $(INSTALL_PATH)/server))
	$(CP) 5point_nister.$(TARGET).$(XSA).xclbin $@
else
	$(warning Could not copy to $(INSTALL_PATH)/server as the folder does not exist)
endif
else
	$(warning Could not find file 5point_nister.$(TARGET).$(XSA).xclbin)
endif

$(XCLBIN_DIR): 
	mkdir -p $@
$(XO_DIR): 
	mkdir -p $@

$(XO_DIR)/5point_nister.xo: $(SRC_DIR)/kernel.cpp | $(XO_DIR)
	$(VPP) $(CLFLAGS) $(LDCLFLAGS) -I $(SRC_DIR) --temp_dir $(XO_DIR) -c -k $(KERNEL_NAME) -o'$@' '$<'
5point_nister.$(TARGET).$(XSA).xclbin: $(KERNEL_XO) | $(XCLBIN_DIR)
	$(VPP) $(CLFLAGS) $(LDCLFLAGS) -I $(SRC_DIR) --temp_dir $(XCLBIN_DIR) --kernel_frequency $(FREQUENCY) -l -o'$@' '$<'

clean:
	-$(RMDIR) _xclbin* _xo* .Xil
	-$(RM) *.log *.jou  *.info *.ltx *.pb *.link_summary

cleanall: clean
	-$(RM) *.xclbin
