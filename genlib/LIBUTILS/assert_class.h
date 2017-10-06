
/******************************************************************************************************************************/
#define HPoffset 8
#define PlacementNewOffset 8

#ifdef _DEBUG
 #define Assert(condition, message)\
{\
  if(!(condition)) std::cerr << (message) << std::endl;\
}
#define Assert2(condition, message, message1)\
{\
  if(!(condition)) std::cerr << (message) << " "<< (message1)<< std::endl;\
}
#define Assert3(condition, message, message1, message2)\
{\
  if(!(condition)) std::cerr << (message) << " " << (message1) << " "<< message2 << std::endl;\
}
#define Assert5(condition, message, message1, message2, message3, message4)\
{\
  if(!(condition)) std::cerr <<(message)<<" "<<(message1)<<" "<<(message2)<<" "<<(message3)<<" "<<(message4)<< std::endl;\
}
 #define _LOG(x) x
 #define CHECK(x) x
#else /* NO_ARG_CHECK */
 #define Assert(condition, message)
 #define Assert2(condition, message, message1)
 #define Assert3(condition, message, message1, message2)
 #define Assert5(condition, message, message1, message2, message3, message4)
 #define _LOG(x)
 #define CHECK(x)
#endif /* NO_ARG_CHECK */
/******************************************************************************************************************************/

