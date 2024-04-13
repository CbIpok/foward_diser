#ifndef TSUNAMI_IOCTL_H
#define TSUNAMI_IOCTL_H

#include <linux/ioctl.h>

#define TSUNAMI_IOCTL_MAGIC 't'

struct tsunamic_dma_xfer {
    unsigned long long src_address;
    unsigned long long dst_address;
    unsigned long long size;
    int dma_handle;
};

#define tsunamic_wait_interrupt _IOR(TSUNAMI_IOCTL_MAGIC, 0, unsigned int)
#define tsunamic_start_dma_tx _IOWR(TSUNAMI_IOCTL_MAGIC, 1, struct tsunamic_dma_xfer)
#define tsunamic_start_dma_rx _IOWR(TSUNAMI_IOCTL_MAGIC, 2, struct tsunamic_dma_xfer)
#define tsunamic_wait_dma _IOW(TSUNAMI_IOCTL_MAGIC, 3, int)

#endif /* TSUNAMI_IOCTL_H */
