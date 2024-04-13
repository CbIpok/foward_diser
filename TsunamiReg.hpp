
#pragma once

#include <cstdint>
#include <string>

namespace forward {
namespace Tsunami {

    class TsunamiControl {
    public:
        static constexpr uint32_t Address() { return 0x00; }

        uint32_t nx : 14; // X dimension
        uint32_t ny : 14; // Y dimension
        uint32_t _unused2928 : 2;
        uint32_t loadDx : 1; // 0 - keep previous x steps
                             // 1 - load new x steps
        uint32_t start : 1; // r - Proccessor ready
                            // w - Start processor

        TsunamiControl()
        {
            *(reinterpret_cast<uint32_t*>(this)) = 0x00000000;
            nx = 0x0;
            ny = 0x0;
            loadDx = 0x0;
            start = 0x0;
        }

        uint32_t pack() const { return *(reinterpret_cast<const uint32_t*>(this)); }

        void unpack(uint32_t value) { *(reinterpret_cast<uint32_t*>(this)) = value; }

        std::string to_string() const
        {
            std::string str = "{";
            str += ".nx = " + std::to_string(nx) + ", ";
            str += ".ny = " + std::to_string(ny) + ", ";
            str += ".loadDx = " + std::to_string(loadDx) + ", ";
            str += ".start = " + std::to_string(start);
            str += "}";
            return str;
        }
    };

    class TsunamiDY {
    public:
        static constexpr uint32_t Address() { return 0x01; }

        uint32_t dy : 32;

        TsunamiDY()
        {
            *(reinterpret_cast<uint32_t*>(this)) = 0x00000000;
            dy = 0x0;
        }

        uint32_t pack() const { return *(reinterpret_cast<const uint32_t*>(this)); }

        void unpack(uint32_t value) { *(reinterpret_cast<uint32_t*>(this)) = value; }

        std::string to_string() const
        {
            std::string str = "{";
            str += ".dy = " + std::to_string(dy);
            str += "}";
            return str;
        }
    };

    class TsunamiDT {
    public:
        static constexpr uint32_t Address() { return 0x02; }

        uint32_t dt : 32;

        TsunamiDT()
        {
            *(reinterpret_cast<uint32_t*>(this)) = 0x00000000;
            dt = 0x0;
        }

        uint32_t pack() const { return *(reinterpret_cast<const uint32_t*>(this)); }

        void unpack(uint32_t value) { *(reinterpret_cast<uint32_t*>(this)) = value; }

        std::string to_string() const
        {
            std::string str = "{";
            str += ".dt = " + std::to_string(dt);
            str += "}";
            return str;
        }
    };

    class TsunamiGround {
    public:
        static constexpr uint32_t Address() { return 0x02; }

        uint32_t grnd : 32;

        TsunamiGround()
        {
            *(reinterpret_cast<uint32_t*>(this)) = 0x00000000;
            grnd = 0x0;
        }

        uint32_t pack() const { return *(reinterpret_cast<const uint32_t*>(this)); }

        void unpack(uint32_t value) { *(reinterpret_cast<uint32_t*>(this)) = value; }

        std::string to_string() const
        {
            std::string str = "{";
            str += ".grnd = " + std::to_string(grnd);
            str += "}";
            return str;
        }
    };

    class IRQEnable {
    public:
        static constexpr uint32_t Address() { return 0x04; }

        uint32_t procIRQ : 1; // 0 - interrupt disabled
                              // 1 - interrupt enabled
        uint32_t _unused3101 : 31;

        IRQEnable()
        {
            *(reinterpret_cast<uint32_t*>(this)) = 0x00000000;
            procIRQ = 0x0;
        }

        uint32_t pack() const { return *(reinterpret_cast<const uint32_t*>(this)); }

        void unpack(uint32_t value) { *(reinterpret_cast<uint32_t*>(this)) = value; }

        std::string to_string() const
        {
            std::string str = "{";
            str += ".procIRQ = " + std::to_string(procIRQ);
            str += "}";
            return str;
        }
    };

    class IRQFlags {
    public:
        static constexpr uint32_t Address() { return 0x05; }

        uint32_t procIRQ : 1;
        uint32_t _unused3101 : 31;

        IRQFlags()
        {
            *(reinterpret_cast<uint32_t*>(this)) = 0x00000000;
            procIRQ = 0x0;
        }

        uint32_t pack() const { return *(reinterpret_cast<const uint32_t*>(this)); }

        void unpack(uint32_t value) { *(reinterpret_cast<uint32_t*>(this)) = value; }

        std::string to_string() const
        {
            std::string str = "{";
            str += ".procIRQ = " + std::to_string(procIRQ);
            str += "}";
            return str;
        }
    };

    class CSR {
    public:
        typedef TsunamiControl TsunamiControlType;
        typedef TsunamiDY TsunamiDYType;
        typedef TsunamiDT TsunamiDTType;
        typedef TsunamiGround TsunamiGroundType;
        typedef IRQEnable IRQEnableType;
        typedef IRQFlags IRQFlagsType;
    };
}
}
